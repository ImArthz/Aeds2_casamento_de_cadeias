#ifndef CASAMENTO_CADEIAS_HPP
#define CASAMENTO_CADEIAS_HPP

#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>

struct ResultadoBusca {
    std::vector<size_t> ocorrencias;
    uint64_t comparacoes_caracteres = 0;
    uint64_t colisoes_espurias = 0;
    double tempo_preprocessamento_ms = 0.0;
    double tempo_busca_ms = 0.0;
    size_t memoria_auxiliar_bytes = 0;
};

struct SeedHit {
    size_t pos_padrao;
    size_t pos_texto;
};

// Algoritmo Ingenuo (Forca Bruta) - Cormen Cap. 32.1
ResultadoBusca busca_ingenua(const std::string& T, const std::string& P);

// Algoritmo de Rabin-Karp - Cormen Cap. 32.2
// d: base do alfabeto (ex: 4 para DNA mapeado ou 256 para ASCII)
// q: numero primo para aritmetica modular
ResultadoBusca rabin_karp_matcher(const std::string& T, const std::string& P, int64_t d = 256, int64_t q = 1000000007LL);

// Funcao de Prefixo (Tabela pi / LPS) do KMP - Cormen Cap. 32.4
std::vector<int> calcular_funcao_prefixo(const std::string& P, uint64_t& comparacoes_pre);

// Algoritmo de Knuth-Morris-Pratt (KMP) - Cormen Cap. 32.4
ResultadoBusca kmp_matcher(const std::string& T, const std::string& P);

// Aplicacao em Bioinformatica: Etapa de Seeding (Semeacao de k-mers) via Rolling Hash
// Utilizado como filtro inicial antes do alinhamento local (estilo BLAST / Seed-and-Extend)
std::vector<SeedHit> seeding_rabin_karp(const std::string& genoma_ref, const std::string& read_query, size_t k, uint64_t& total_comparacoes);

#endif
