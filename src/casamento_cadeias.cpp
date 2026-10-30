#include "casamento_cadeias.hpp"
#include <chrono>
#include <unordered_map>

using Relogio = std::chrono::high_resolution_clock;

ResultadoBusca busca_ingenua(const std::string& T, const std::string& P) {
    ResultadoBusca res;
    res.memoria_auxiliar_bytes = sizeof(size_t) * 2;

    size_t n = T.size();
    size_t m = P.size();
    if (m == 0 || m > n) {
        return res;
    }

    auto inicio = Relogio::now();

    for (size_t s = 0; s <= n - m; ++s) {
        size_t j = 0;
        while (j < m) {
            res.comparacoes_caracteres++;
            if (T[s + j] != P[j]) {
                break;
            }
            j++;
        }
        if (j == m) {
            res.ocorrencias.push_back(s);
        }
    }

    auto fim = Relogio::now();
    res.tempo_busca_ms = std::chrono::duration<double, std::milli>(fim - inicio).count();
    return res;
}

ResultadoBusca rabin_karp_matcher(const std::string& T, const std::string& P, int64_t d, int64_t q) {
    ResultadoBusca res;
    res.memoria_auxiliar_bytes = sizeof(int64_t) * 5;

    size_t n = T.size();
    size_t m = P.size();
    if (m == 0 || m > n) {
        return res;
    }

    auto inicio_pre = Relogio::now();

    // h = d^(m-1) mod q
    int64_t h = 1;
    for (size_t i = 0; i < m - 1; ++i) {
        h = (h * d) % q;
    }

    int64_t p = 0;
    int64_t t = 0;

    // Calculo inicial dos hashes p e t_0 pela Regra de Horner
    for (size_t i = 0; i < m; ++i) {
        p = (d * p + static_cast<unsigned char>(P[i])) % q;
        t = (d * t + static_cast<unsigned char>(T[i])) % q;
    }

    auto fim_pre = Relogio::now();
    res.tempo_preprocessamento_ms = std::chrono::duration<double, std::milli>(fim_pre - inicio_pre).count();

    auto inicio_busca = Relogio::now();

    for (size_t s = 0; s <= n - m; ++s) {
        if (p == t) {
            // Verificacao caractere a caractere para descartar colisao espuria
            size_t j = 0;
            while (j < m) {
                res.comparacoes_caracteres++;
                if (T[s + j] != P[j]) {
                    break;
                }
                j++;
            }
            if (j == m) {
                res.ocorrencias.push_back(s);
            } else {
                res.colisoes_espurias++;
            }
        }

        // Atualiza o rolling hash para t_{s+1} em O(1)
        if (s < n - m) {
            int64_t char_sai = static_cast<unsigned char>(T[s]);
            int64_t char_entra = static_cast<unsigned char>(T[s + m]);

            t = (d * (t - char_sai * h) + char_entra) % q;

            // Em C++, o operador % pode retornar negativo se o dividendo for negativo
            if (t < 0) {
                t += q;
            }
        }
    }

    auto fim_busca = Relogio::now();
    res.tempo_busca_ms = std::chrono::duration<double, std::milli>(fim_busca - inicio_busca).count();
    return res;
}

std::vector<int> calcular_funcao_prefixo(const std::string& P, uint64_t& comparacoes_pre) {
    size_t m = P.size();
    std::vector<int> pi(m, 0);
    int k = 0;

    for (size_t q = 1; q < m; ++q) {
        while (k > 0) {
            comparacoes_pre++;
            if (P[k] == P[q]) {
                break;
            }
            k = pi[k - 1];
        }
        comparacoes_pre++;
        if (P[k] == P[q]) {
            k++;
        }
        pi[q] = k;
    }
    return pi;
}

ResultadoBusca kmp_matcher(const std::string& T, const std::string& P) {
    ResultadoBusca res;
    size_t n = T.size();
    size_t m = P.size();
    if (m == 0 || m > n) {
        return res;
    }

    auto inicio_pre = Relogio::now();
    uint64_t comp_pre = 0;
    std::vector<int> pi = calcular_funcao_prefixo(P, comp_pre);
    auto fim_pre = Relogio::now();

    res.tempo_preprocessamento_ms = std::chrono::duration<double, std::milli>(fim_pre - inicio_pre).count();
    res.memoria_auxiliar_bytes = pi.capacity() * sizeof(int);
    res.comparacoes_caracteres += comp_pre;

    auto inicio_busca = Relogio::now();
    int q = 0; // numero de caracteres casados

    for (size_t i = 0; i < n; ++i) {
        while (q > 0) {
            res.comparacoes_caracteres++;
            if (P[q] == T[i]) {
                break;
            }
            q = pi[q - 1];
        }
        res.comparacoes_caracteres++;
        if (P[q] == T[i]) {
            q++;
        }
        if (static_cast<size_t>(q) == m) {
            res.ocorrencias.push_back(i + 1 - m);
            q = pi[q - 1];
        }
    }

    auto fim_busca = Relogio::now();
    res.tempo_busca_ms = std::chrono::duration<double, std::milli>(fim_busca - inicio_busca).count();
    return res;
}

std::vector<SeedHit> seeding_rabin_karp(const std::string& genoma_ref, const std::string& read_query, size_t k, uint64_t& total_comparacoes) {
    std::vector<SeedHit> seeds;
    total_comparacoes = 0;

    if (k == 0 || read_query.size() < k || genoma_ref.size() < k) {
        return seeds;
    }

    const int64_t d = 4;
    const int64_t q = 1000000007LL;

    auto mapear_base = [](char c) -> int64_t {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 0;
        }
    };

    int64_t h = 1;
    for (size_t i = 0; i < k - 1; ++i) {
        h = (h * d) % q;
    }

    // Indexa os k-mers da read usando rolling hash
    std::unordered_map<int64_t, std::vector<size_t>> indice_kmer;
    int64_t hash_read = 0;
    for (size_t i = 0; i < k; ++i) {
        hash_read = (d * hash_read + mapear_base(read_query[i])) % q;
    }
    indice_kmer[hash_read].push_back(0);

    for (size_t s = 0; s < read_query.size() - k; ++s) {
        int64_t sai = mapear_base(read_query[s]);
        int64_t entra = mapear_base(read_query[s + k]);
        hash_read = (d * (hash_read - sai * h) + entra) % q;
        if (hash_read < 0) hash_read += q;
        indice_kmer[hash_read].push_back(s + 1);
    }

    // Varre o genoma de referencia em O(n) buscando seeds exatas
    int64_t hash_ref = 0;
    for (size_t i = 0; i < k; ++i) {
        hash_ref = (d * hash_ref + mapear_base(genoma_ref[i])) % q;
    }

    for (size_t s = 0; s <= genoma_ref.size() - k; ++s) {
        auto it = indice_kmer.find(hash_ref);
        if (it != indice_kmer.end()) {
            for (size_t pos_r : it->second) {
                bool igual = true;
                for (size_t j = 0; j < k; ++j) {
                    total_comparacoes++;
                    if (genoma_ref[s + j] != read_query[pos_r + j]) {
                        igual = false;
                        break;
                    }
                }
                if (igual) {
                    seeds.push_back({pos_r, s});
                }
            }
        }

        if (s < genoma_ref.size() - k) {
            int64_t sai = mapear_base(genoma_ref[s]);
            int64_t entra = mapear_base(genoma_ref[s + k]);
            hash_ref = (d * (hash_ref - sai * h) + entra) % q;
            if (hash_ref < 0) hash_ref += q;
        }
    }

    return seeds;
}
