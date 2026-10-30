import os
import csv
import matplotlib.pyplot as plt
import numpy as np

# Configuracao visual padrao artigo academico (IEEE / SBC)
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "axes.grid": True,
    "grid.alpha": 0.35,
    "grid.linestyle": "--"
})

DIR_BENCH = "benchmarks"
DIR_OUT = "docs/relatorio_latex"
os.makedirs(DIR_OUT, exist_ok=True)


def ler_csv(caminho):
    with open(caminho, mode="r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def plotar_escala_n():
    dados = ler_csv(os.path.join(DIR_BENCH, "resultados_escala_n.csv"))
    algoritmos = ["Ingenuo", "Rabin-Karp", "KMP"]
    cores = {"Ingenuo": "#d95f02", "Rabin-Karp": "#7570b3", "KMP": "#1b9e77"}
    marcadores = {"Ingenuo": "o", "Rabin-Karp": "s", "KMP": "^"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    for alg in algoritmos:
        linhas = [d for d in dados if d["Algoritmo"] == alg]
        x = [int(d["TamanhoN"]) / 1e6 for d in linhas]
        y_tempo = [float(d["TempoTotalMs"]) for d in linhas]
        y_comp = [int(d["Comparacoes"]) / 1e6 for d in linhas]

        ax1.plot(x, y_tempo, label=alg, color=cores[alg], marker=marcadores[alg], linewidth=1.8)
        ax2.plot(x, y_comp, label=alg, color=cores[alg], marker=marcadores[alg], linewidth=1.8)

    ax1.set_title("(a) Tempo Total de Execucao vs. Tamanho do Genoma (N)")
    ax1.set_xlabel("Tamanho da Sequencia T (Milhoes de pares de base - Mpb)")
    ax1.set_ylabel("Tempo de Execucao (ms)")
    ax1.legend()

    ax2.set_title("(b) Comparacoes de Caracteres vs. Tamanho do Genoma (N)")
    ax2.set_xlabel("Tamanho da Sequencia T (Milhoes de pares de base - Mpb)")
    ax2.set_ylabel("Comparacoes Diretas de Bases (Milhoes)")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(DIR_OUT, "tempo_escala_n.png"), dpi=300)
    plt.close()


def plotar_pior_caso():
    dados = ler_csv(os.path.join(DIR_BENCH, "resultados_pior_caso.csv"))
    algoritmos = ["Ingenuo", "Rabin-Karp", "KMP"]
    cores = {"Ingenuo": "#d95f02", "Rabin-Karp": "#7570b3", "KMP": "#1b9e77"}
    marcadores = {"Ingenuo": "o", "Rabin-Karp": "s", "KMP": "^"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    for alg in algoritmos:
        linhas = [d for d in dados if d["Algoritmo"] == alg]
        x = [int(d["TamanhoM"]) for d in linhas]
        y_tempo = [float(d["TempoTotalMs"]) for d in linhas]
        y_comp = [int(d["Comparacoes"]) for d in linhas]

        ax1.plot(x, y_tempo, label=alg, color=cores[alg], marker=marcadores[alg], linewidth=1.8)
        ax2.plot(x, y_comp, label=alg, color=cores[alg], marker=marcadores[alg], linewidth=1.8)

    ax1.set_title("(a) Tempo no Pior Caso (DNA Poli-A, N = 100.000 pb)")
    ax1.set_xlabel("Tamanho do Padrao Repetitivo M (pb)")
    ax1.set_ylabel("Tempo de Execucao (ms) - Escala Log")
    ax1.set_yscale("log")
    ax1.legend()

    ax2.set_title("(b) Comparacoes de Caracteres no Pior Caso (Escala Log)")
    ax2.set_xlabel("Tamanho do Padrao Repetitivo M (pb)")
    ax2.set_ylabel("Total de Comparacoes (Escala Log)")
    ax2.set_yscale("log")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(DIR_OUT, "pior_caso_repetitivo.png"), dpi=300)
    plt.close()


def plotar_colisoes_rk():
    dados = ler_csv(os.path.join(DIR_BENCH, "resultados_colisoes_rk.csv"))
    primos = [f"q={d['PrimoQ']}" for d in dados]
    colisoes = [int(d["ColisoesEspurias"]) for d in dados]
    tempos = [float(d["TempoBuscaMs"]) for d in dados]

    fig, ax1 = plt.subplots(figsize=(8.5, 4.4))

    x = np.arange(len(primos))
    largura = 0.42

    barras = ax1.bar(x - largura / 2, colisoes, width=largura, color="#7570b3", label="Colisoes Espurias (Falsos Hits de Hash)")
    ax1.set_xlabel("Modulo Primo Escolhido (q) no Rabin-Karp (N = 500.000 pb, M = 32 pb)")
    ax1.set_ylabel("Numero de Colisoes Espurias (Escala Log)", color="#4a4480")
    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels(primos, rotation=20)

    for b, c in zip(barras, colisoes):
        altura = max(c, 1)
        ax1.annotate(f"{c}", xy=(b.get_x() + b.get_width() / 2, altura),
                     xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=8)

    ax2 = ax1.twinx()
    ax2.plot(x + largura / 2, tempos, color="#d95f02", marker="o", linewidth=2.0, label="Tempo de Busca (ms)")
    ax2.set_ylabel("Tempo de Busca (ms)", color="#d95f02")
    ax2.grid(False)

    linhas1, labels1 = ax1.get_legend_handles_labels()
    linhas2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(linhas1 + linhas2, labels1 + labels2, loc="upper right")

    plt.title("Impacto da Escolha do Primo q nas Colisoes Espurias e Tempo do Rabin-Karp")
    plt.tight_layout()
    plt.savefig(os.path.join(DIR_OUT, "colisoes_rabin_karp.png"), dpi=300)
    plt.close()


def plotar_genoma_real():
    dados = ler_csv(os.path.join(DIR_BENCH, "resultados_genoma_real.csv"))
    padroes = []
    for d in dados:
        if d["Padrao"] not in padroes:
            padroes.append(d["Padrao"])

    algs = ["Ingenuo", "Rabin-Karp_q1e9", "Rabin-Karp_q101", "KMP"]
    nomes_legenda = {
        "Ingenuo": "Ingenuo (Forca Bruta)",
        "Rabin-Karp_q1e9": "Rabin-Karp (q = 10^9+7)",
        "Rabin-Karp_q101": "Rabin-Karp (q = 101)",
        "KMP": "Knuth-Morris-Pratt (KMP)"
    }
    cores = {
        "Ingenuo": "#d95f02",
        "Rabin-Karp_q1e9": "#7570b3",
        "Rabin-Karp_q101": "#e7298a",
        "KMP": "#1b9e77"
    }

    x = np.arange(len(padroes))
    largura = 0.2

    fig, ax = plt.subplots(figsize=(9.5, 4.6))

    for idx, alg in enumerate(algs):
        valores = []
        for p in padroes:
            item = next(d for d in dados if d["Padrao"] == p and d["Algoritmo"] == alg)
            valores.append(max(int(item["Comparacoes"]), 1))
        desloc = (idx - 1.5) * largura
        ax.bar(x + desloc, valores, width=largura, label=nomes_legenda[alg], color=cores[alg])

    ax.set_title("Comparacoes Diretas de Bases no Genoma Real do SARS-CoV-2 (NC_045512.2 - 29.903 pb)")
    ax.set_xlabel("Padrao Biologico Buscado")
    ax.set_ylabel("Numero de Comparacoes de Caracteres (Escala Log)")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(padroes)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(DIR_OUT, "comparacoes_genoma_sars.png"), dpi=300)
    plt.close()


if __name__ == "__main__":
    plotar_escala_n()
    plotar_pior_caso()
    plotar_colisoes_rk()
    plotar_genoma_real()
    print("Graficos gerados em docs/relatorio_latex/")
