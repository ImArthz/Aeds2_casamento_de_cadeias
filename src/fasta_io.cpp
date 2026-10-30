#include "fasta_io.hpp"
#include <fstream>
#include <stdexcept>
#include <random>
#include <cctype>

SequenciaFasta ler_arquivo_fasta(const std::string& caminho) {
    std::ifstream arq(caminho);
    if (!arq.is_open()) {
        throw std::runtime_error("Nao foi possivel abrir o arquivo FASTA: " + caminho);
    }

    SequenciaFasta resultado;
    std::string linha;
    resultado.sequencia.reserve(32768);

    while (std::getline(arq, linha)) {
        if (linha.empty()) continue;
        if (linha[0] == '>') {
            if (resultado.cabecalho.empty()) {
                resultado.cabecalho = linha.substr(1);
            }
            continue;
        }
        for (char c : linha) {
            if (!std::isspace(static_cast<unsigned char>(c))) {
                resultado.sequencia.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
            }
        }
    }

    return resultado;
}

std::string gerar_dna_aleatorio(size_t tamanho, unsigned int semente) {
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    std::mt19937 gerador(semente);
    std::uniform_int_distribution<int> dist(0, 3);

    std::string dna;
    dna.resize(tamanho);
    for (size_t i = 0; i < tamanho; ++i) {
        dna[i] = bases[dist(gerador)];
    }
    return dna;
}

std::string gerar_dna_repetitivo(size_t tamanho_texto, size_t tamanho_padrao, std::string& padrao_saida) {
    // Texto composto por 'A's com 'C' no final para forcar (n-m+1)*m comparacoes no algoritmo ingenuo
    std::string texto(tamanho_texto, 'A');
    if (tamanho_texto > 0) {
        texto[tamanho_texto - 1] = 'C';
    }

    padrao_saida.assign(tamanho_padrao, 'A');
    if (tamanho_padrao > 0) {
        padrao_saida[tamanho_padrao - 1] = 'C';
    }

    return texto;
}
