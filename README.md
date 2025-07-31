# Metodos-Simplex-Matlab
Este repositório contém uma coleção de scripts em MATLAB que implementam diferentes variações do método Simplex para resolver problemas de programação linear (PL): Uma Fase, Duas Fases, Artificial Única e Big-M. Cada script aborda um cenário específico e inclui um exemplo prático para demonstrar sua aplicação.

## Implementações do Método Simplex em MATLAB
Este repositório contém uma coleção de scripts em MATLAB que implementam diferentes variações do método Simplex (Uma Fase, Duas Fases, Artificial Única e Big-M) para resolver problemas de programação linear (PL). Cada script aborda um cenário específico e inclui um exemplo prático para demonstrar sua aplicação.

<img width="409" height="249" alt="image" src="https://github.com/user-attachments/assets/4aaa1a53-07f4-4707-94d2-6086639f2e5a" />

Métodos Implementados:
1.  **Simplex Uma Fase (`uma_fase.m`)**
2.  **Simplex Duas Fases (`duas_fases.m`)**
3.  **Método Big M (`bigm.m`)**
4.  **Método da Artificial Única (`artificial_unica.m`)**

---

## O Método Simplex de Uma Fase
O Simplex de Uma Fase é a aplicação mais direta do algoritmo Simplex. Ele é utilizado em problemas de programação linear onde uma solução básica factível inicial é facilmente identificada. Geralmente, isso ocorre em problemas cujas restrições são todas do tipo "menor ou igual" ($\le$), pois as variáveis de folga introduzidas para converter as inequações em igualdades formam uma base inicial viável (matriz identidade).

### Implementação: `uma_fase.m`
Este script implementa o método Simplex tabular para resolver um PPL que já possui uma base inicial factível. O código inicializa a tabela Simplex e itera, selecionando a variável para entrar na base e a variável para sair da base a cada passo, até que a condição de otimalidade seja alcançada (quando não há custos reduzidos que possam melhorar a função objetivo).

#### Exemplo Resolvido no Código
O script está configurado para resolver o seguinte problema de minimização, já no formato padrão (com as variáveis de folga $x_3, x_4, x_5, x_6$ formando a base inicial):

$$
\begin{aligned}
\text{Minimizar } Z = 5x_1 &+ 4x_2 \\
\text{Sujeito a:} \\
6x_1 + 4x_2 + x_3 &= 24 \\
x_1 + 2x_2 + x_4 &= 6 \\
-x_1 + x_2 + x_5 &= 1 \\
x_2 + x_6 &= 2 \\
x_i &\ge 0, \quad \forall i \in \{1, ..., 6\}
\end{aligned}
$$

---

## O Método Simplex de Duas Fases
O Método de Duas Fases é uma abordagem robusta utilizada quando uma solução básica factível não é óbvia, o que é comum em problemas com restrições de igualdade ($=$) ou "maior ou igual" ($\ge$).

- **Fase 1:** O objetivo é encontrar uma solução básica factível para o problema original. Para isso, introduz-se variáveis artificiais nas restrições que não possuem uma variável básica inicial. Em seguida, resolve-se um novo PPL cujo objetivo é minimizar a soma dessas variáveis artificiais. Se ao final da Fase 1 o valor ótimo for zero, significa que uma solução factível foi encontrada e todas as variáveis artificiais são nulas.
- **Fase 2:** A tabela final da Fase 1 (com as colunas das variáveis artificiais removidas e a função objetivo original restaurada) é usada como a tabela inicial para otimizar o problema original.

### Implementação: `duas_fases.m`
Este script implementa o algoritmo Simplex de Duas Fases. Ele automatiza a criação da função objetivo artificial para a Fase 1, resolve-a para encontrar uma base factível e, em caso de sucesso, transiciona para a Fase 2 para otimizar a função objetivo original do problema.

#### Exemplo Resolvido no Código
O script resolve o seguinte problema de minimização. As variáveis $x_6$ e $x_7$ são as variáveis artificiais adicionadas para iniciar a Fase 1, enquanto $x_5$ é uma variável de folga.

$$
\begin{aligned}
\text{Minimizar } Z = -x_1 &+ 2x_2 \\
\text{Sujeito a:} \\
x_1 + x_2 - x_3 &= 2 \\
-x_1 + x_2 - x_4 &= 1 \\
x_2 + x_5 &= 3 \\
x_i &\ge 0, \quad \forall i \in \{1, ..., 5\}
\end{aligned}
$$

---

## O Método da Variável Artificial Única
Esta é uma técnica engenhosa para lidar com problemas que, ao serem colocados na forma padrão, apresentam valores negativos no vetor de recursos (RHS), tornando a base inicial infactível. O método consiste em introduzir uma única variável artificial ($x_a$) em todas as restrições com RHS negativo. Em seguida, um pivoteamento especial é realizado na coluna dessa variável artificial e na linha correspondente à restrição mais infactível (aquela com o RHS mais negativo) para tornar a tabela factível. A partir daí, o processo segue de forma similar ao Método de Duas Fases.

### Implementação: `artificial_unica.m`
O código `artificial_unica.m` aplica esta técnica. Primeiramente, ele gera uma tabela Simplex inicial que é infactível devido a valores negativos no RHS. Em seguida, realiza um pivoteamento estratégico na coluna da variável artificial para tornar a solução básica factível. Após essa etapa, o algoritmo prossegue com a Fase 1 para eliminar a variável artificial da base e então para a Fase 2, otimizando o problema original.

#### Exemplo Resolvido no Código
O script foi projetado para um problema cujas restrições originais levam a um RHS negativo na forma padrão. A variável $x_5$ é a única variável artificial utilizada para restaurar a factibilidade.

$$
\begin{aligned}
\text{Minimizar } Z = -2x_1 &- 3x_2 \\
\text{Sujeito a:} \\
-x_1 - x_2 + x_3 &= -3 \\
2x_1 - x_2 + x_4 &= -2 \\
x_i &\ge 0, \quad \forall i \in \{1, ..., 4\}
\end{aligned}
$$

---

## O Método Big M
O Método Big M é uma alternativa ao Método de Duas Fases para lidar com problemas que necessitam de variáveis artificiais. Em vez de separar o processo em duas fases, o Big M incorpora as variáveis artificiais diretamente na função objetivo original. A cada variável artificial é associado um custo de penalidade muito grande, representado por "M". Em um problema de minimização, as variáveis artificiais entram com um custo de $+M$; em maximização, com um custo de $-M$. Essa penalidade força o algoritmo a zerar essas variáveis na solução ótima, se uma solução factível existir.

### Implementação: `bigm.m`
Este script resolve um PPL usando o método Big M. Ele utiliza a capacidade simbólica do MATLAB para manipular a variável `M`. O código constrói a tabela Simplex inicial com os custos `M` associados às variáveis artificiais na função objetivo e prossegue com as iterações do Simplex até encontrar a solução ótima.

#### Exemplo Resolvido no Código
O problema exemplo requer variáveis artificiais ($x_8, x_9$) para as restrições de "maior ou igual" e uma variável de folga ($x_5$) para a restrição de "menor ou igual". As variáveis $x_6$ e $x_7$ são as variáveis de excesso.

$$
\begin{aligned}
\text{Minimizar } Z = 2x_1 &- 2x_2 - x_3 - x_4 \\
\text{Sujeito a:} \\
x_1 + 2x_2 + 2x_3 + x_4 + x_5 &= 2 \\
x_1 - 2x_2 + x_3 + 2x_4 - x_6 + x_8 &= 3 \\
2x_1 - x_2 + 3x_3 - x_7 + x_9 &= 2 \\
x_i &\ge 0, \text{ para todas as variáveis}
\end{aligned}
$$
