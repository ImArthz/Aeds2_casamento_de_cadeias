#include "casamento_cadeias.hpp"
#include "fasta_io.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>

void imprimir_linha_tabela(const std::string& padrao_nome, const std::string& alg, const ResultadoBusca& r) {
    std::cout << std::left << std::setw(24) << padrao_nome
              << std::setw(18) << alg
              << std::right << std::setw(10) << r.ocorrencias.size()
              << std::setw(16) << r.comparacoes_caracteres
              << std::setw(14) << r.colisoes_espurias
              << std::setw(14) << std::fixed << std::setprecision(4) << (r.tempo_preprocessamento_ms + r.tempo_busca_ms)
              << "\n";
}

void executar_estudo_caso_ncbi() {
    std::cout << "========================================================================================\n";
    std::cout << "1. ESTUDO DE CASO EM BIOINFORMATICA: GENOMAS REAIS (NCBI NC_045512.2 e NM_000546.6)\n";
    std::cout << "========================================================================================\n";

    SequenciaFasta sars = ler_arquivo_fasta("data/sars_cov2_NC_045512.fasta");
    SequenciaFasta tp53 = ler_arquivo_fasta("data/tp53_NM_000546.fasta");

    std::cout << "Genoma carregado [1]: " << sars.cabecalho << " (" << sars.sequencia.size() << " pb)\n";
    std::cout << "Genoma carregado [2]: " << tp53.cabecalho << " (" << tp53.sequencia.size() << " pb)\n\n";

    struct CasoBio {
        std::string nome;
        std::string padrao;
    };

    std::vector<CasoBio> casos = {
        {"Motif TATA-box", "TATAAA"},
        {"Primer OMS RdRp-F2", "GTGAAATGGTCATGTGTGGCGG"},
        {"Sitio Furina Spike", "CCTCGGCGGGCACGTAGT"},
        {"Cauda Poli-A (25pb)", "AAAAAAAAAAAAAAAAAAAAAAAAA"}
    };

    std::ofstream arq_bio("benchmarks/resultados_genoma_real.csv");
    arq_bio << "Dataset,Padrao,TamanhoM,Algoritmo,Ocorrencias,Comparacoes,ColisoesEspurias,TempoPreMs,TempoBuscaMs,TempoTotalMs\n";

    std::cout << std::left << std::setw(24) << "Padrao Biologico"
              << std::setw(18) << "Algoritmo"
              << std::right << std::setw(10) << "Hits"
              << std::setw(16) << "Comparacoes"
              << std::setw(14) << "Colisoes"
              << std::setw(14) << "Tempo (ms)" << "\n";
    std::cout << "----------------------------------------------------------------------------------------\n";

    for (const auto& c : casos) {
        ResultadoBusca r_ing = busca_ingenua(sars.sequencia, c.padrao);
        ResultadoBusca r_rk_bom = rabin_karp_matcher(sars.sequencia, c.padrao, 4, 1000000007LL);
        ResultadoBusca r_rk_ruim = rabin_karp_matcher(sars.sequencia, c.padrao, 4, 101LL);
        ResultadoBusca r_kmp = kmp_matcher(sars.sequencia, c.padrao);

        assert(r_ing.ocorrencias == r_rk_bom.ocorrencias);
        assert(r_ing.ocorrencias == r_kmp.ocorrencias);

        imprimir_linha_tabela(c.nome, "Ingenuo", r_ing);
        imprimir_linha_tabela(c.nome, "Rabin-Karp (q=10^9+7)", r_rk_bom);
        imprimir_linha_tabela(c.nome, "Rabin-Karp (q=101)", r_rk_ruim);
        imprimir_linha_tabela(c.nome, "KMP", r_kmp);
        std::cout << "----------------------------------------------------------------------------------------\n";

        auto gravar = [&](const std::string& alg, const ResultadoBusca& r) {
            arq_bio << "SARS-CoV-2," << c.nome << "," << c.padrao.size() << "," << alg << ","
                    << r.ocorrencias.size() << "," << r.comparacoes_caracteres << ","
                    << r.colisoes_espurias << "," << r.tempo_preprocessamento_ms << ","
                    << r.tempo_busca_ms << "," << (r.tempo_preprocessamento_ms + r.tempo_busca_ms) << "\n";
        };
        gravar("Ingenuo", r_ing);
        gravar("Rabin-Karp_q1e9", r_rk_bom);
        gravar("Rabin-Karp_q101", r_rk_ruim);
        gravar("KMP", r_kmp);
    }

    // Demonstracao de Seeding (Seed-and-Extend) com leitura contendo mutacoes pontuais (SNPs)
    size_t pos_real = 21562; // Inicio do gene Spike (S) no SARS-CoV-2
    std::string read_simulada = sars.sequencia.substr(pos_real, 120);
    // Introduz 2 mutacoes pontuais (SNPs) que impediriam casamento exato da read inteira
    read_simulada[35] = (read_simulada[35] == 'A') ? 'G' : 'A';
    read_simulada[88] = (read_simulada[88] == 'T') ? 'C' : 'T';

    ResultadoBusca busca_direta = kmp_matcher(sars.sequencia, read_simulada);
    uint64_t comp_seed = 0;
    std::vector<SeedHit> seeds = seeding_rabin_karp(sars.sequencia, read_simulada, 15, comp_seed);

    std::cout << "\n[Aplicacao em Alinhamento Genetico - Seeding de k-mers (k=15)]\n";
    std::cout << "Read de 120 pb extraida da posicao " << pos_real << " (Gene Spike) com 2 mutacoes SNP.\n";
    std::cout << "-> Casamento exato da read inteira (KMP): " << busca_direta.ocorrencias.size() << " ocorrencias.\n";
    std::cout << "-> Seeding via Rolling Hash (Rabin-Karp k=15): " << seeds.size()
              << " ancoras exatas encontradas (diagonal principal em pos=" << (seeds.empty() ? 0 : seeds[0].pos_texto) << ").\n\n";
}

void executar_benchmarks_assintoticos() {
    std::cout << "========================================================================================\n";
    std::cout << "2. EXECUTANDO BENCHMARKS ASSINTOTICOS (ESCALA N, PIOR CASO E COLISOES RABIN-KARP)\n";
    std::cout << "========================================================================================\n";

    // Experimento A: Variando tamanho do texto N (DNA aleatorio, M = 64)
    {
        std::ofstream arq("benchmarks/resultados_escala_n.csv");
        arq << "TamanhoN,TamanhoM,Algoritmo,Comparacoes,ColisoesEspurias,TempoPreMs,TempoBuscaMs,TempoTotalMs,MemoriaBytes\n";

        std::vector<size_t> tamanhos_n = {10000, 50000, 100000, 250000, 500000, 1000000, 2500000, 5000000};
        size_t m = 64;

        for (size_t n : tamanhos_n) {
            std::string T = gerar_dna_aleatorio(n, 1337);
            std::string P = T.substr(n - m, m);

            ResultadoBusca r_ing = busca_ingenua(T, P);
            ResultadoBusca r_rk = rabin_karp_matcher(T, P, 4, 1000000007LL);
            ResultadoBusca r_kmp = kmp_matcher(T, P);

            assert(r_ing.ocorrencias == r_kmp.ocorrencias);

            auto salvar = [&](const std::string& nome, const ResultadoBusca& r) {
                arq << n << "," << m << "," << nome << ","
                    << r.comparacoes_caracteres << "," << r.colisoes_espurias << ","
                    << r.tempo_preprocessamento_ms << "," << r.tempo_busca_ms << ","
                    << (r.tempo_preprocessamento_ms + r.tempo_busca_ms) << ","
                    << r.memoria_auxiliar_bytes << "\n";
            };
            salvar("Ingenuo", r_ing);
            salvar("Rabin-Karp", r_rk);
            salvar("KMP", r_kmp);
            std::cout << "[Escala N] Concluido N = " << n << " pb\n";
        }
    }

    // Experimento B: Pior Caso Patologico em DNA Repetitivo (Poli-A: AAAAA...AC)
    {
        std::ofstream arq("benchmarks/resultados_pior_caso.csv");
        arq << "TamanhoN,TamanhoM,Algoritmo,Comparacoes,ColisoesEspurias,TempoPreMs,TempoBuscaMs,TempoTotalMs\n";

        size_t n = 100000;
        std::vector<size_t> tamanhos_m = {50, 100, 250, 500, 1000, 2000};

        for (size_t m : tamanhos_m) {
            std::string P;
            std::string T = gerar_dna_repetitivo(n, m, P);

            ResultadoBusca r_ing = busca_ingenua(T, P);
            ResultadoBusca r_rk = rabin_karp_matcher(T, P, 4, 1000000007LL);
            ResultadoBusca r_kmp = kmp_matcher(T, P);

            assert(r_ing.ocorrencias == r_kmp.ocorrencias);

            auto salvar = [&](const std::string& nome, const ResultadoBusca& r) {
                arq << n << "," << m << "," << nome << ","
                    << r.comparacoes_caracteres << "," << r.colisoes_espurias << ","
                    << r.tempo_preprocessamento_ms << "," << r.tempo_busca_ms << ","
                    << (r.tempo_preprocessamento_ms + r.tempo_busca_ms) << "\n";
            };
            salvar("Ingenuo", r_ing);
            salvar("Rabin-Karp", r_rk);
            salvar("KMP", r_kmp);
            std::cout << "[Pior Caso Repetitivo] Concluido N = " << n << ", M = " << m << " pb\n";
        }
    }

    // Experimento C: Impacto do Modulo Primo q nas Colisoes Espurias do Rabin-Karp
    {
        std::ofstream arq("benchmarks/resultados_colisoes_rk.csv");
        arq << "PrimoQ,TamanhoN,TamanhoM,OcorrenciasReais,ColisoesEspurias,Comparacoes,TempoBuscaMs\n";

        size_t n = 500000;
        size_t m = 32;
        std::string T = gerar_dna_aleatorio(n, 2026);
        std::string P = T.substr(n / 2, m);

        std::vector<int64_t> primos = {13, 31, 101, 503, 1009, 10007, 100003, 1000000007LL};
        for (int64_t q : primos) {
            ResultadoBusca r_rk = rabin_karp_matcher(T, P, 4, q);
            arq << q << "," << n << "," << m << ","
                << r_rk.ocorrencias.size() << "," << r_rk.colisoes_espurias << ","
                << r_rk.comparacoes_caracteres << "," << r_rk.tempo_busca_ms << "\n";
            std::cout << "[Colisoes RK] q = " << std::setw(10) << q
                      << " -> Colisoes espurias: " << r_rk.colisoes_espurias << "\n";
        }
    }
}

int main() {
    executar_estudo_caso_ncbi();
    executar_benchmarks_assintoticos();
    std::cout << "\nTodos os experimentos foram concluidos e exportados para a pasta benchmarks/.\n";
    return 0;
}
