# Roteiro da Aula Técnica (Pitch de 40 Minutos)
**Disciplina:** Algoritmos e Estruturas de Dados II / Projeto e Análise de Algoritmos  
**Professor:** Michel Pires da Silva (DECOM-DV – CEFET-MG)  
**Tema 5:** Casamento de Cadeias (*String Matching* - Cap. 32 do Cormen): Algoritmos de Rabin-Karp, Knuth-Morris-Pratt (KMP) e Implicações em Bioinformática

---

## Divisão da Equipe (5 Integrantes — 8 Minutos por Integrante)

- **Integrante 1 (00:00 – 08:00):** `c[NOME DO INTEGRANTE 1]` — *Problema Motivador em Genômica, Formalização (Cormen Cap. 32.1) e Por Que a Força Bruta Falha no DNA.*
- **Integrante 2 (08:00 – 16:00):** `[NOME DO INTEGRANTE 2]` — *Algoritmo de Rabin-Karp (Cap. 32.2): Regra de Horner, Rolling Hash em $O(1)$ e Colisões Espúrias.*
- **Integrante 3 (16:00 – 24:00):** `[NOME DO INTEGRANTE 3]` — *Algoritmo Knuth-Morris-Pratt (Cap. 32.4): Propriedade Sufixo-Prefixo, Tabela $\pi$ e Demonstração Visual no Simulador.*
- **Integrante 4 (24:00 – 32:00):** `[NOME DO INTEGRANTE 4]` — *Análise Assintótica Rigorosa ($O, \Omega, \Theta$), Prova Amortizada do KMP e Aplicação no Alinhamento Genético (Seeding).*
- **Integrante 5 (32:00 – 40:00):** `[NOME DO INTEGRANTE 5]` — *Engenharia da Implementação em C++17, Gargalos de Hardware/Cache nos Testes Reais (NCBI SARS-CoV-2) e Debate com a Turma.*

---

## Cronograma Detalhado Passo a Passo

### Bloco 1 (00:00 – 08:00) — `[NOME DO INTEGRANTE 1]`
**Objetivo:** Instigar a turma com um problema real de bioinformática e mostrar no quadro onde o algoritmo ingênuo desperdiça trabalho.

1. **Pergunta Motivadora para a Sala (00:00 – 02:30):**
   - *"Imagine que acabamos de sequenciar uma amostra clínica e temos um arquivo FASTA com o genoma humano ($n \approx 3,2 \times 10^9$ pares de base) ou milhares de leituras virais de SARS-CoV-2. Precisamos localizar um primer de PCR de $m = 22$ bases ou uma sequência promotora TATA-box. Se testarmos posição por posição voltando do zero a cada erro, quantas operações faremos?"*
2. **Formalização no Quadro (02:30 – 05:00):**
   - Definir texto $T[0..n-1]$ e padrão $P[0..m-1]$ sobre o alfabeto genético $\Sigma = \{\text{A}, \text{C}, \text{G}, \text{T}\}$ com $|\Sigma| = 4$.
   - Definir deslocamento válido (*valid shift*) $s \in [0, n-m]$ tal que $T[s..s+m-1] = P[0..m-1]$.
3. **A Armadilha do Alfabeto Pequeno ($|\Sigma| = 4$) na Busca Ingênua (05:00 – 08:00):**
   - Escrever no quadro uma região de repetição em tandem real (cauda poli-A do SARS-CoV-2: `T = AAAAAAAAAAAAAAC`, `P = AAAAAC`).
   - Mostrar que a probabilidade de duas bases aleatórias coincidirem no DNA é $1/4 = 25\%$ (muito maior que no texto em português, que é $1/26 \approx 3,8\%$).
   - Em regiões repetitivas, o algoritmo ingênuo compara `AAAAA`, erra no `C`, volta o ponteiro do texto para `s+1` e **re-compara todos os `A`s que já sabia que estavam lá**, caindo em $\Theta((n-m+1)m)$.

---

### Bloco 2 (08:00 – 16:00) — `[NOME DO INTEGRANTE 2]`
**Objetivo:** Construir a ideia algébrica do Rabin-Karp passo a passo: transformar subcadeias de DNA em números inteiros.

1. **De Letras para Polinômios — Regra de Horner (08:00 – 10:30):**
   - Mapear $\text{A}=0, \text{C}=1, \text{G}=2, \text{T}=3$ (base $d = 4$).
   - Mostrar como uma janela de tamanho $m$ vira um número: $p = P[0]\cdot d^{m-1} + P[1]\cdot d^{m-2} + \dots + P[m-1]$.
   - Com a Regra de Horner, calculamos $p$ e a primeira janela $t_0$ em $\Theta(m)$ multiplicações.
2. **O Salto de Engenharia: *Rolling Hash* em $\Theta(1)$ (10:30 – 13:30):**
   - Questionar a turma: *"Se recalcularmos o polinômio do zero a cada passo, continuamos em $\Theta(nm)$. Como aproveitar $t_s$ para achar $t_{s+1}$ em tempo constante?"*
   - Desenhar no quadro a janela deslizando:
     1. Subtrai o caractere que sai pela esquerda: $t_s - T[s]\cdot d^{m-1}$;
     2. Multiplica tudo pela base $d$ (deslocamento à esquerda);
     3. Soma o novo caractere que entra pela direita: $+ T[s+m]$.
   - Aplicar o módulo primo $q$ para evitar *overflow* de inteiros de 64 bits:
     $$t_{s+1} = \left(d \cdot (t_s - T[s]\cdot h) + T[s+m]\right) \bmod q, \quad \text{onde } h = d^{m-1} \bmod q$$
3. **Colisões Espúrias (*Spurious Hits*) (13:30 – 16:00):**
   - Explicar que $t_s \equiv p \pmod q$ **não garante** que as strings são iguais (princípio da casa dos pombos).
   - Quando os hashes batem, fazemos a verificação caractere a caractere. Se $q$ for grande ($q = 10^9 + 7$), a taxa de colisão é $\approx 1/q$, mantendo o caso médio em $O(n+m)$.

---

### Bloco 3 (16:00 – 24:00) — `[NOME DO INTEGRANTE 3]`
**Objetivo:** Apresentar o KMP como a solução determinística que **nunca retrocede o ponteiro $i$ do texto** e demonstrar ao vivo no Simulador Visual.

1. **A Intuição de Knuth, Morris e Pratt (16:00 – 18:30):**
   - Voltar ao exemplo onde casamos 5 caracteres de `P` e erramos no 6º.
   - Se já lemos $T[s..s+q-1]$ e sabemos que esse trecho é idêntico a $P[0..q-1]$, nós já conhecemos exatamente quais letras estão no texto!
   - Só precisamos saber: **qual é o maior prefixo próprio de $P[0..q-1]$ que também é sufixo de $P[0..q-1]$?**
2. **Construção da Tabela de Prefixos $\pi[0..m-1]$ (18:30 – 21:30):**
   - Fazer no quadro o exemplo `P = ACAGACACAG`:
     - Mostrar como $\pi[q]$ guarda o tamanho do maior prefixo de $P$ que é sufixo próprio de $P_q$.
     - Explicar por que quando ocorre falha na posição $q$, fazemos `q = pi[q - 1]` (em indexação base-0) sem tocar no índice $i$ do texto.
3. **Demonstração Ao Vivo no Simulador Visual (21:30 – 24:00):**
   - Abrir o arquivo `docs/seminario/simulador_interativo.html` no projetor.
   - Executar o **Preset 1 (Poli-A)** e o **Preset 2 (Sobreposição KMP)** passo a passo, mostrando a diferença entre as comparações da Busca Ingênua e do KMP.

---

### Bloco 4 (24:00 – 32:00) — `[NOME DO INTEGRANTE 4]`
**Objetivo:** Formalizar a análise assintótica rigorosa ($O, \Omega, \Theta$) e explicar como Rabin-Karp e KMP entram no Alinhamento Genético real.

1. **Debate de Análise Assintótica com a Turma (24:00 – 28:00):**
   - **Pergunta Provocativa para a Turma:** *"Olhem para o código do KMP: temos um laço `for` de $1$ até $n$ e, dentro dele, um laço `while (q > 0 && P[q] != T[i])`. Dois laços aninhados não deveriam dar $O(nm)$? Por que o KMP é $\Theta(n)$?"*
   - **Explicação Amortizada no Quadro (Método Potencial/Agregado do Cormen):**
     - A variável de estado $q$ começa em $0$.
     - A cada iteração do `for`, $q$ pode aumentar **no máximo em 1 unidade** (`q++` ocorre no máximo $n$ vezes em todo o programa).
     - Cada iteração do `while` interno substitui $q$ por $\pi[q-1] < q$, ou seja, **diminui $q$ em pelo menos 1 unidade**.
     - Como $q$ nunca fica negativo ($q \ge 0$), o `while` não pode diminuir $q$ mais vezes do que ele foi incrementado! Logo, o corpo do `while` executa no máximo $n$ vezes somando **todas** as iterações do `for`. Custo total: $\Theta(m)$ no pré-processamento + $\Theta(n)$ na busca = $\Theta(n+m)$ no pior caso!
2. **Conexão com Bioinformática: Casamento Exato vs. Alinhamento Genético (28:00 – 32:00):**
   - Diferenciar **Casamento Exato** ($\Theta(n+m)$) de **Alinhamento Global/Local com Mutações** (Needleman-Wunsch / Smith-Waterman, que exige Programação Dinâmica em $\Theta(nm)$).
   - Mostrar que alinhar uma leitura (*read*) contra um genoma inteiro via Programação Dinâmica é computacionalmente inviável.
   - Explicar o paradigma **Seed-and-Extend (usado no BLAST e alinhadores modernos)**:
     1. **Seeding (Rabin-Karp / $k$-mers):** Usa *rolling hash* em $\Theta(n)$ para encontrar pequenas sementes exatas (*seeds* de tamanho $k=11$ a $15$) mesmo quando a sequência tem mutações SNP em outras posições;
     2. **Extending:** Roda programação dinâmica apenas ao redor das âncoras encontradas pelo casamento exato!

---

### Bloco 5 (32:00 – 40:00) — `[NOME DO INTEGRANTE 5]`
**Objetivo:** Apresentar os resultados reais em C++17 com o genoma do SARS-CoV-2 (`NC_045512.2`), discutir os gargalos de hardware/linguagem e fechar com o Cheat Sheet.

1. **Resultados Práticos no Genoma do SARS-CoV-2 e Escala Progressiva (32:00 – 35:30):**
   - Mostrar os gráficos gerados a partir do código em C++17 (`tempo_escala_n.png`, `pior_caso_repetitivo.png` e `colisoes_rabin_karp.png`).
   - Destacar o teste de *Seeding*: quando injetamos 2 mutações SNP em uma leitura de 120 pb do gene *Spike*, a busca exata da leitura inteira retorna 0 ocorrências, mas o *seeding* por Rabin-Karp ($k=15$) recupera 91 âncoras exatas apontando diretamente para o locus correto (`pos = 21.562`).
2. **Gargalos Críticos de Hardware e Linguagem Descobertos na Prática (35:30 – 38:30):**
   - **Por que o Rabin-Karp tem constante oculta mais alta que o KMP em C++ para um único padrão?**
     - Na ALU do processador x86_64, a instrução de divisão inteira/resto (`idiv` gerada pelo operador `% q`) custa de **20 a 40 ciclos de clock**, enquanto a comparação de bytes na cache L1 (`cmp`) do KMP custa **1 ciclo de clock**.
   - **Armadilhas de Implementação em C++:**
     - O operador `%` em C/C++ preserva o sinal do dividendo: `(t - T[s]*h) % q` pode resultar em número negativo, exigindo a correção `if (t < 0) t += q;`.
     - Necessidade de inteiros de 64 bits (`int64_t`) para evitar *integer overflow* na multiplicação $d \cdot t_s$ antes do módulo.
3. **Fechamento e Perguntas Finais (38:30 – 40:00):**
   - Apresentar o **Cheat Sheet** (quando usar Força Bruta, quando usar Rabin-Karp para múltiplos padrões/*seeding*, e quando usar KMP para garantia linear $\Theta(n+m)$ em alfabetos pequenos como o DNA).
