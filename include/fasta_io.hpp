#ifndef FASTA_IO_HPP
#define FASTA_IO_HPP

#include <string>
#include <vector>
#include <cstddef>

struct SequenciaFasta {
    std::string cabecalho;
    std::string sequencia;
};

// Le um arquivo FASTA (.fasta / .fna) concatenando as linhas da sequencia genetica
SequenciaFasta ler_arquivo_fasta(const std::string& caminho);

// Gera sequencia de DNA pseudo-aleatoria uniforme sobre o alfabeto {A, C, G, T}
std::string gerar_dna_aleatorio(size_t tamanho, unsigned int semente = 42);

// Gera sequencia de DNA repetitiva (ex: repeticoes em tandem poli-A "AAAAAAA...AC")
// Usada para demonstrar a degradacao quadratica O(nm) da busca ingenua no alfabeto de DNA
std::string gerar_dna_repetitivo(size_t tamanho_texto, size_t tamanho_padrao, std::string& padrao_saida);

#endif
