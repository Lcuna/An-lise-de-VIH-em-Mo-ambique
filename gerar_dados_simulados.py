"""
gerar_dados_simulados.py
========================
Gera dados FICTÍCIOS de vigilância do VIH em Moçambique (2019–2025),
para construir e testar o script de análise epidemiológica.

ATENÇÃO: todos os valores são simulados. Servem apenas para treino e
não devem ser interpretados como dados reais do país.

Ficheiros criados em dados/brutos/:
  1. casos_vih.csv        – lista nominal (1 linha = 1 caso notificado)
  2. ptv_mensal.csv       – relatório mensal agregado da PTV por província
  3. populacao.csv        – denominadores: população por província/ano/sexo/idade
  4. parametros_verdadeiros.json        – os "valores reais" usados na simulação
                                          (servem para validar a análise mais tarde)
  5. registo_problemas_introduzidos.json – erros de qualidade inseridos de propósito
                                          (o Passo 2, limpeza, deve encontrá-los)

Uso:  python gerar_dados_simulados.py
"""

import json
import sys
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):      # acentos na consola do Windows
    sys.stdout.reconfigure(encoding="utf-8")

SEMENTE = 2026                      # garante que os dados são sempre os mesmos
rng = np.random.default_rng(SEMENTE)

PASTA_PROJETO = Path(__file__).resolve().parent
PASTA_BRUTOS = PASTA_PROJETO / "dados" / "brutos"

ANOS = list(range(2019, 2026))
ANO_BASE = 2019

# ---------------------------------------------------------------------------
# 1. PARÂMETROS DA SIMULAÇÃO ("a verdade" que a análise deverá recuperar)
# ---------------------------------------------------------------------------

# Província: (população total em 2019, risco relativo de diagnóstico de VIH)
PROVINCIAS = {
    "Niassa":           (1_900_000, 0.55),
    "Cabo Delgado":     (2_400_000, 0.60),
    "Nampula":          (6_000_000, 0.45),
    "Zambézia":         (5_300_000, 0.95),
    "Tete":             (2_800_000, 0.70),
    "Manica":           (2_000_000, 1.10),
    "Sofala":           (2_300_000, 1.15),
    "Inhambane":        (1_500_000, 1.05),
    "Gaza":             (1_450_000, 1.80),
    "Maputo Província": (2_000_000, 1.70),
    "Maputo Cidade":    (1_150_000, 1.40),
}
CRESCIMENTO_POP = 0.028             # crescimento populacional anual
PROP_MASCULINA = 0.48

# Grupo etário: (idade mín, idade máx, proporção da população,
#                RR feminino vs. referência, RR masculino vs. feminino)
# Referência = mulheres 25–29 anos. O grupo 0–4 é gerado a partir da PTV.
GRUPOS_ETARIOS = {
    "0-4":   (0, 4,   0.160, None, None),
    "5-9":   (5, 9,   0.145, 0.04, 1.00),
    "10-14": (10, 14, 0.125, 0.06, 1.00),
    "15-19": (15, 19, 0.110, 0.45, 0.25),
    "20-24": (20, 24, 0.095, 0.95, 0.35),
    "25-29": (25, 29, 0.080, 1.00, 0.55),
    "30-34": (30, 34, 0.065, 0.90, 0.75),
    "35-39": (35, 39, 0.053, 0.75, 0.90),
    "40-44": (40, 44, 0.042, 0.60, 1.00),
    "45-49": (45, 49, 0.033, 0.45, 1.00),
    "50+":   (50, 79, 0.092, 0.20, 1.05),
}

TAXA_BASE_100MIL = 350              # mulheres 25–29, RR provincial = 1, ano 2019
TENDENCIA_ANUAL = -0.05             # log(RR) por ano → ≈ -4,9 % ao ano
EFEITO_COVID = {2020: 0.78, 2021: 0.92}   # quebra na testagem


def peso_mensal(ano: int) -> np.ndarray:
    """Distribuição dos casos pelos 12 meses (sazonalidade + COVID em 2020)."""
    meses = np.arange(1, 13)
    pesos = 1 + 0.15 * np.cos(2 * np.pi * (meses - 7) / 12)   # pico em julho
    if ano == 2020:
        pesos[3:7] *= 0.55          # abril–julho 2020: confinamento
    return pesos / pesos.sum()


# PTV (prevenção da transmissão vertical)
TAXA_NATALIDADE = 0.038
COBERTURA_CPN = 0.95
PREV_GESTANTES_BASE = 0.09          # multiplicada pelo RR provincial
TRANSMISSAO_COM_TARV = 0.02
TRANSMISSAO_SEM_TARV = 0.18


def interpolar(ano: int, inicio: float, fim: float) -> float:
    """Valor linear entre 2019 (inicio) e 2025 (fim)."""
    return inicio + (fim - inicio) * (ano - ANOS[0]) / (ANOS[-1] - ANOS[0])


# ---------------------------------------------------------------------------
# 2. POPULAÇÃO (denominadores)
# ---------------------------------------------------------------------------
# Estrutura etária diferente por província: >1 = população mais jovem (mais crianças,
# menos idosos), <1 = mais envelhecida/urbana. É isto que torna a padronização útil.
FATOR_JUVENTUDE = {
    "Niassa": 1.12, "Cabo Delgado": 1.08, "Nampula": 1.08, "Zambézia": 1.06,
    "Tete": 1.08, "Manica": 1.00, "Sofala": 0.95, "Inhambane": 0.92,
    "Gaza": 0.90, "Maputo Província": 0.80, "Maputo Cidade": 0.65,
}


def estrutura_etaria(prov: str) -> dict:
    f = FATOR_JUVENTUDE[prov]
    pesos = {}
    for grupo, (idade_min, _, prop, _, _) in GRUPOS_ETARIOS.items():
        if idade_min < 15:
            pesos[grupo] = prop * f
        elif idade_min >= 50:
            pesos[grupo] = prop / f
        else:
            pesos[grupo] = prop
    total = sum(pesos.values())
    return {g: p / total for g, p in pesos.items()}


def gerar_populacao() -> pd.DataFrame:
    linhas = []
    for prov, (pop2019, _) in PROVINCIAS.items():
        estrutura = estrutura_etaria(prov)
        for ano in ANOS:
            pop_total = pop2019 * (1 + CRESCIMENTO_POP) ** (ano - ANO_BASE)
            for grupo, prop in estrutura.items():
                for sexo, p_sexo in (("M", PROP_MASCULINA), ("F", 1 - PROP_MASCULINA)):
                    linhas.append({
                        "provincia": prov, "ano": ano, "sexo": sexo,
                        "grupo_etario": grupo,
                        "populacao": int(round(pop_total * prop * p_sexo)),
                    })
    return pd.DataFrame(linhas)


# ---------------------------------------------------------------------------
# 3. RELATÓRIO MENSAL DA PTV
# ---------------------------------------------------------------------------
def gerar_ptv(pop: pd.DataFrame) -> pd.DataFrame:
    pop_prov = pop.groupby(["provincia", "ano"])["populacao"].sum()
    linhas = []
    for prov, (_, rr_prov) in PROVINCIAS.items():
        ruido_prov = rng.normal(0, 0.02)          # cada província tem o seu nível
        for ano in ANOS:
            prev = PREV_GESTANTES_BASE * rr_prov * (0.97 ** (ano - ANO_BASE))
            cob_teste = interpolar(ano, 0.93, 0.98)
            cob_tarv = np.clip(interpolar(ano, 0.85, 0.97) + ruido_prov, 0, 0.995)
            cob_pcr = interpolar(ano, 0.65, 0.85)
            fator_covid = EFEITO_COVID.get(ano, 1.0)
            for mes in range(1, 13):
                esperadas = pop_prov[(prov, ano)] * TAXA_NATALIDADE / 12 * COBERTURA_CPN
                cpn = rng.poisson(esperadas * (0.9 if ano == 2020 and 4 <= mes <= 7 else 1))
                testadas = rng.binomial(cpn, cob_teste * (0.95 if fator_covid < 1 else 1))
                positivas = rng.binomial(testadas, prev)
                em_tarv = rng.binomial(positivas, cob_tarv)
                expostas = rng.binomial(positivas, 0.95)
                pcr = rng.binomial(expostas, cob_pcr)
                prop_tarv = em_tarv / positivas if positivas else 1
                p_transm = (prop_tarv * TRANSMISSAO_COM_TARV
                            + (1 - prop_tarv) * TRANSMISSAO_SEM_TARV)
                pcr_pos = rng.binomial(pcr, p_transm)
                linhas.append({
                    "provincia": prov, "ano": ano, "mes": mes,
                    "gestantes_cpn": cpn,
                    "gestantes_testadas": testadas,
                    "gestantes_positivas": positivas,
                    "gestantes_em_tarv": em_tarv,
                    "criancas_expostas": expostas,
                    "criancas_pcr_2meses": pcr,
                    "criancas_pcr_positivo": pcr_pos,
                })
    return pd.DataFrame(linhas)


# ---------------------------------------------------------------------------
# 4. LISTA NOMINAL DE CASOS
# ---------------------------------------------------------------------------
PONTOS_ENTRADA = {
    "F_adulta":  (["ATS comunitário", "Consulta externa", "CPN", "Maternidade", "TB", "Outros"],
                  [0.28, 0.30, 0.20, 0.05, 0.05, 0.12]),
    "F_idosa":   (["ATS comunitário", "Consulta externa", "TB", "Outros"],
                  [0.35, 0.45, 0.08, 0.12]),
    "M":         (["ATS comunitário", "Consulta externa", "TB", "Outros"],
                  [0.40, 0.40, 0.12, 0.08]),
    "crianca":   (["Consulta externa", "ATS comunitário", "TB", "Internamento pediátrico"],
                  [0.50, 0.20, 0.10, 0.20]),
}


def datas_aleatorias(ano: int, n: int, mes: int | None = None) -> pd.DatetimeIndex:
    """Sorteia n datas dentro de um ano (ou de um mês)."""
    if mes is None:
        meses = rng.choice(np.arange(1, 13), size=n, p=peso_mensal(ano))
    else:
        meses = np.full(n, mes)
    dias = rng.integers(1, 29, size=n)          # 1–28: válido em todos os meses
    return pd.to_datetime(pd.DataFrame({"year": ano, "month": meses, "day": dias}))


def gerar_casos(pop: pd.DataFrame, ptv: pd.DataFrame) -> pd.DataFrame:
    blocos = []

    # 4a. Casos ≥ 5 anos: modelo de Poisson (taxa × população)
    for linha in pop.itertuples(index=False):
        idade_min, idade_max, _, rr_f, rr_m = GRUPOS_ETARIOS[linha.grupo_etario]
        if rr_f is None:
            continue
        rr_idade = rr_f * (rr_m if linha.sexo == "M" else 1)
        rr_prov = PROVINCIAS[linha.provincia][1]
        rr_ano = np.exp(TENDENCIA_ANUAL * (linha.ano - ANO_BASE)) * EFEITO_COVID.get(linha.ano, 1)
        esperado = linha.populacao * TAXA_BASE_100MIL / 1e5 * rr_idade * rr_prov * rr_ano
        n = rng.poisson(esperado)
        if n == 0:
            continue
        blocos.append(pd.DataFrame({
            "provincia": linha.provincia, "sexo": linha.sexo,
            "idade": rng.integers(idade_min, idade_max + 1, size=n),
            "data_diagnostico": datas_aleatorias(linha.ano, n),
        }))

    # 4b. Crianças 0–4 anos: vêm da PTV
    for linha in ptv.itertuples(index=False):
        # bebés com PCR positivo (idade 0), diagnosticados na CCR
        n0 = linha.criancas_pcr_positivo
        if n0:
            blocos.append(pd.DataFrame({
                "provincia": linha.provincia,
                "sexo": rng.choice(["M", "F"], size=n0),
                "idade": 0,
                "data_diagnostico": datas_aleatorias(linha.ano, n0, linha.mes),
                "ponto_entrada": "CCR (PCR criança exposta)",
            }))
        # crianças infetadas que escaparam ao PCR → diagnóstico tardio (1–4 anos)
        n_tardio = rng.poisson(linha.criancas_pcr_positivo * 0.5)
        if n_tardio:
            blocos.append(pd.DataFrame({
                "provincia": linha.provincia,
                "sexo": rng.choice(["M", "F"], size=n_tardio),
                "idade": rng.integers(1, 5, size=n_tardio),
                "data_diagnostico": datas_aleatorias(linha.ano, n_tardio),
                "ponto_entrada": rng.choice(["Consulta externa", "Internamento pediátrico"],
                                            size=n_tardio, p=[0.6, 0.4]),
            }))

    casos = pd.concat(blocos, ignore_index=True)
    n = len(casos)

    # Ponto de entrada (para quem ainda não tem)
    sem_pe = casos["ponto_entrada"].isna() if "ponto_entrada" in casos else pd.Series(True, index=casos.index)
    for chave, (opcoes, probs) in PONTOS_ENTRADA.items():
        if chave == "crianca":
            mascara = sem_pe & (casos["idade"] < 15)
        elif chave == "M":
            mascara = sem_pe & (casos["idade"] >= 15) & (casos["sexo"] == "M")
        elif chave == "F_adulta":
            mascara = sem_pe & casos["idade"].between(15, 49) & (casos["sexo"] == "F")
        else:
            mascara = sem_pe & (casos["idade"] >= 50) & (casos["sexo"] == "F")
        casos.loc[mascara, "ponto_entrada"] = rng.choice(opcoes, size=mascara.sum(), p=probs)

    # Gestante
    casos["gestante"] = "Não"
    casos.loc[casos["ponto_entrada"].isin(["CPN", "Maternidade"]), "gestante"] = "Sim"
    outras_mulheres = (casos["sexo"] == "F") & casos["idade"].between(15, 49) & (casos["gestante"] == "Não")
    casos.loc[outras_mulheres & (rng.random(n) < 0.05), "gestante"] = "Sim"
    casos.loc[casos["sexo"] == "M", "gestante"] = "Não aplicável"
    casos.loc[casos["idade"] < 10, "gestante"] = "Não aplicável"

    # CD4 no diagnóstico (células/mm³): mais baixo em homens e em idades mais velhas
    mediana = np.where(casos["sexo"] == "M", 300, 380) * np.where(casos["idade"] >= 40, 0.85, 1)
    casos["cd4_diagnostico"] = np.round(np.exp(rng.normal(np.log(mediana), 0.6))).astype(float)
    casos.loc[casos["idade"] < 5, "cd4_diagnostico"] = np.nan        # em bebés usa-se % CD4
    casos.loc[rng.random(n) < 0.35, "cd4_diagnostico"] = np.nan

    # Data de notificação: 0–30 dias depois do diagnóstico
    casos["data_notificacao"] = casos["data_diagnostico"] + pd.to_timedelta(
        rng.integers(0, 31, size=n), unit="D")

    # Início de TARV ("testar e tratar"): a cobertura melhora com os anos
    ano_diag = casos["data_diagnostico"].dt.year
    p_inicia = 0.85 + 0.10 * (ano_diag - 2019) / 6
    inicia = rng.random(n) < p_inicia
    atraso = np.round(rng.exponential(7, size=n))
    casos["data_inicio_tarv"] = casos["data_diagnostico"] + pd.to_timedelta(atraso, unit="D")
    casos.loc[~inicia, "data_inicio_tarv"] = pd.NaT

    # Ordenar, criar ID, baralhar
    casos = casos.sort_values("data_diagnostico").reset_index(drop=True)
    casos.insert(0, "id_caso", [f"MZ-VIH-{i:06d}" for i in range(1, n + 1)])
    return casos


# ---------------------------------------------------------------------------
# 5. INTRODUZIR PROBLEMAS DE QUALIDADE (como nos dados reais!)
# ---------------------------------------------------------------------------
def sem_acentos(texto: str) -> str:
    return "".join(c for c in unicodedata.normalize("NFD", texto) if unicodedata.category(c) != "Mn")


VARIANTES_PROVINCIA = {
    "Maputo Cidade": ["Cidade de Maputo", "MAPUTO CIDADE", "maputo cidade ", "C. Maputo"],
    "Maputo Província": ["Maputo Provincia", "Prov. Maputo", "MAPUTO PROVÍNCIA"],
    "Cabo Delgado": ["C. Delgado", "cabo delgado", "Cabo-Delgado"],
}


def variante_provincia(prov: str) -> str:
    opcoes = VARIANTES_PROVINCIA.get(prov, []) + [prov.upper(), prov.lower(), f" {prov} ", sem_acentos(prov)]
    return str(rng.choice(opcoes))


def sujar_casos(casos: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    df = casos.copy()
    n = len(df)
    reg = {}

    def sortear(p):
        return rng.random(n) < p

    # Datas como texto, para podermos misturar formatos
    for col in ["data_diagnostico", "data_notificacao", "data_inicio_tarv"]:
        df[col] = df[col].dt.strftime("%Y-%m-%d").astype(object)

    m = sortear(0.08)
    df.loc[m, "provincia"] = [variante_provincia(p) for p in df.loc[m, "provincia"]]
    reg["provincia_grafia_inconsistente"] = int(m.sum())

    mapa_m = ["M", "Masculino", "m", "masc"]
    mapa_f = ["F", "Feminino", "f", "fem"]
    m = sortear(0.10)
    df.loc[m & (df["sexo"] == "M"), "sexo"] = rng.choice(mapa_m, size=(m & (df["sexo"] == "M")).sum())
    df.loc[m & (df["sexo"] == "F"), "sexo"] = rng.choice(mapa_f, size=(m & (df["sexo"] == "F")).sum())
    reg["sexo_codificacao_variada"] = int(m.sum())
    m = sortear(0.02)
    df.loc[m, "sexo"] = rng.choice(["", "Desconhecido", None], size=m.sum())
    reg["sexo_em_falta"] = int(m.sum())

    df["idade"] = df["idade"].astype(float)
    m = sortear(0.02)
    df.loc[m, "idade"] = np.nan
    reg["idade_em_falta"] = int(m.sum())
    m = sortear(0.003)
    df.loc[m, "idade"] = rng.choice([-1, 150, 999], size=m.sum())
    reg["idade_impossivel"] = int(m.sum())

    m = sortear(0.20) & df["data_diagnostico"].notna()
    df.loc[m, "data_diagnostico"] = pd.to_datetime(df.loc[m, "data_diagnostico"]).dt.strftime("%d/%m/%Y")
    reg["data_diagnostico_formato_dd_mm_aaaa"] = int(m.sum())
    m = sortear(0.01)
    df.loc[m, "data_diagnostico"] = None
    reg["data_diagnostico_em_falta"] = int(m.sum())

    m = sortear(0.01)
    diag = pd.to_datetime(casos.loc[m, "data_diagnostico"])
    df.loc[m, "data_notificacao"] = (diag - pd.to_timedelta(rng.integers(5, 60, size=m.sum()), unit="D")).dt.strftime("%Y-%m-%d")
    reg["notificacao_antes_do_diagnostico"] = int(m.sum())

    m = sortear(0.003) & df["data_inicio_tarv"].notna()
    df.loc[m, "data_inicio_tarv"] = (pd.to_datetime(casos.loc[m, "data_diagnostico"])
                                     - pd.to_timedelta(rng.integers(10, 200, size=m.sum()), unit="D")).dt.strftime("%Y-%m-%d")
    reg["tarv_antes_do_diagnostico"] = int(m.sum())

    m = sortear(0.005) & (casos["sexo"] == "M")
    df.loc[m, "gestante"] = "Sim"
    reg["homem_registado_como_gestante"] = int(m.sum())

    m = sortear(0.002)
    df.loc[m, "cd4_diagnostico"] = rng.choice([-5, 9999, 15000], size=m.sum())
    reg["cd4_impossivel"] = int(m.sum())

    # Duplicados exatos (mesma ficha introduzida 2 vezes)
    dup = df.sample(frac=0.015, random_state=SEMENTE)
    reg["duplicados_exatos"] = len(dup)
    # Re-notificações (mesmo ID, outra data de notificação)
    renot = df.sample(frac=0.005, random_state=SEMENTE + 1).copy()
    renot["data_notificacao"] = (pd.to_datetime(renot["data_notificacao"], errors="coerce")
                                 + pd.Timedelta(days=45)).dt.strftime("%Y-%m-%d")
    reg["renotificacoes_mesmo_id"] = len(renot)

    df = pd.concat([df, dup, renot], ignore_index=True)
    df = df.sample(frac=1, random_state=SEMENTE).reset_index(drop=True)   # baralhar
    return df, reg


def sujar_ptv(ptv: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    df = ptv.copy()
    reg = {}
    n = len(df)

    falta = rng.random(n) < 0.03
    reg["relatorios_mensais_em_falta"] = int(falta.sum())
    df = df.loc[~falta].reset_index(drop=True)
    n = len(df)

    m = rng.random(n) < 0.05
    df.loc[m, "provincia"] = [variante_provincia(p) for p in df.loc[m, "provincia"]]
    reg["provincia_grafia_inconsistente"] = int(m.sum())

    m = rng.random(n) < 0.01
    df.loc[m, "gestantes_testadas"] = (df.loc[m, "gestantes_cpn"] * 1.1).round().astype(int)
    reg["testadas_maior_que_cpn"] = int(m.sum())

    df = df.astype({c: float for c in df.columns if c.startswith(("gestantes", "criancas"))})
    m = rng.random(n) < 0.01
    cols = rng.choice([c for c in df.columns if c.startswith(("gestantes", "criancas"))], size=m.sum())
    for i, c in zip(df.index[m], cols):
        df.loc[i, c] = np.nan
    reg["valores_em_falta"] = int(m.sum())

    dup = df.sample(n=6, random_state=SEMENTE)
    reg["linhas_duplicadas"] = len(dup)
    df = pd.concat([df, dup], ignore_index=True).sample(frac=1, random_state=SEMENTE).reset_index(drop=True)
    return df, reg


# ---------------------------------------------------------------------------
# 6. EXECUTAR
# ---------------------------------------------------------------------------
def main():
    PASTA_BRUTOS.mkdir(parents=True, exist_ok=True)

    pop = gerar_populacao()
    ptv = gerar_ptv(pop)
    casos = gerar_casos(pop, ptv)

    casos_sujos, reg_casos = sujar_casos(casos)
    ptv_sujo, reg_ptv = sujar_ptv(ptv)

    # utf-8-sig → o Excel no Windows abre os acentos corretamente
    pop.to_csv(PASTA_BRUTOS / "populacao.csv", index=False, encoding="utf-8-sig")
    ptv_sujo.to_csv(PASTA_BRUTOS / "ptv_mensal.csv", index=False, encoding="utf-8-sig")
    casos_sujos.to_csv(PASTA_BRUTOS / "casos_vih.csv", index=False, encoding="utf-8-sig")

    parametros = {
        "aviso": "Parâmetros usados para SIMULAR os dados (fictícios).",
        "semente": SEMENTE,
        "taxa_base_por_100mil": TAXA_BASE_100MIL,
        "referencia": "mulheres 25-29 anos, RR provincial = 1, ano 2019",
        "rr_provincial": {p: v[1] for p, v in PROVINCIAS.items()},
        "rr_idade_mulheres": {g: v[3] for g, v in GRUPOS_ETARIOS.items() if v[3] is not None},
        "rr_homens_vs_mulheres": {g: v[4] for g, v in GRUPOS_ETARIOS.items() if v[4] is not None},
        "tendencia_anual_log_rr": TENDENCIA_ANUAL,
        "variacao_anual_percent": round((np.exp(TENDENCIA_ANUAL) - 1) * 100, 2),
        "efeito_covid_rr": EFEITO_COVID,
        "transmissao_vertical": {"mae_em_tarv": TRANSMISSAO_COM_TARV,
                                 "mae_sem_tarv": TRANSMISSAO_SEM_TARV},
        "n_casos_unicos": len(casos),
    }
    with open(PASTA_BRUTOS / "parametros_verdadeiros.json", "w", encoding="utf-8") as f:
        json.dump(parametros, f, ensure_ascii=False, indent=2)
    with open(PASTA_BRUTOS / "registo_problemas_introduzidos.json", "w", encoding="utf-8") as f:
        json.dump({"casos_vih": reg_casos, "ptv_mensal": reg_ptv}, f, ensure_ascii=False, indent=2)

    print("Dados simulados criados em:", PASTA_BRUTOS)
    print(f"  casos_vih.csv   : {len(casos_sujos):>7,} linhas ({len(casos):,} casos únicos)")
    print(f"  ptv_mensal.csv  : {len(ptv_sujo):>7,} linhas")
    print(f"  populacao.csv   : {len(pop):>7,} linhas")


if __name__ == "__main__":
    main()
