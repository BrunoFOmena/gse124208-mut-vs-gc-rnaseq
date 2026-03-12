## Projeto

Modelo exploratório de separação transcriptômica entre células ERCC6 mutantes (síndrome de Cockayne) e gene‑corrected (GC) no conjunto público GSE124208, utilizando redução de dimensionalidade por PCA e regressão logística em um pipeline supervisionado com validação cruzada.

### Introdução

O conjunto GSE124208 reúne dados de RNA‑seq de células derivadas de um paciente com síndrome de Cockayne (mutação em ERCC6) e de suas contrapartes gene‑corrigidas (GC). As amostras abrangem células iPSC, MSC e NSC expostas ou não a radiação UV e estresse replicativo.  
Esses dados permitem investigar se o estado mutante versus gene‑corrigido se reflete em perfis transcriptômicos diferenciáveis. O objetivo deste projeto não é construir um classificador diagnóstico clínico, mas sim avaliar, de forma exploratória, até que ponto o transcriptoma separa mut de GC em um espaço de baixa dimensão sob um modelo supervisionado simples.

### Objetivo do projeto

Construir e avaliar um **modelo exploratório de separação transcriptômica mut vs GC** a partir das contagens NCBI‑generated do GSE124208, integradas ao series matrix, utilizando:
- filtragem de genes por expressão mínima;
- transformação `log2(CPM+1)` e PCA;
- regressão logística com validação cruzada estratificada;
com foco na **caracterização da separação de classes** e não na performance diagnóstica generalizável.

### Metodologia

- **Dados de entrada**
  - Contagens brutas NCBI (`GSE124208_raw_counts_GRCh38.p13_NCBI.tsv`): genes × amostras.
  - Metadados (`GSE124208_series_matrix.txt`): rótulos mut/GC, tipo celular (iPSC, MSC, NSC) e condição (Ctrl, UV, RS).
  - Anotações de genes (`Human.GRCh38.p13.annot.tsv`) para interpretação opcional.

- **Pré‑processamento**
  - Construção de um `DataFrame` de metadados com `sample_id`, `status` (mut/GC), `cell_type` e `condition`.  
  - Alinhamento exato entre colunas de contagens e amostras dos metadados, com conferência registrada em `TABELAS/conferencia_amostras_contagens.csv`.
  - Filtro **não supervisionado** de baixa expressão: remoção de genes com poucas contagens em poucas amostras (`TABELAS/filtro_genes_resumo.csv` e `genes_filtrados.csv`).
  - Transformação de contagens para `CPM` e, em seguida, `X_log = log2(CPM+1)`. Estatísticas globais antes/depois são registradas em `TABELAS/normalizacao_resumo.csv`.

- **PCA exploratório**
  - PCA global em `X_log` para fins descritivos, gerando:
    - `TABELAS/X_pca_exploratorio.csv` (scores por amostra, com metadados);
    - `TABELAS/pca_loadings_exploratorio.csv` (loadings por gene e componente).
  - Visualização de PC1 vs PC2 colorida por status mut/GC (`FIGURAS/pca_mut_vs_gc.png`), para inspeção da separação bruta no espaço de componentes.

- **Modelo supervisionado e validação**
  - Pipeline `StandardScaler → PCA (5 componentes) → LogisticRegression (L2, solver liblinear, max_iter=1000)`.
  - Validação cruzada estratificada em 5 folds (`StratifiedKFold`, shuffle, `random_state=42`), usando apenas `X_log` e o vetor binário Y (mut=1, GC=0).
  - Em cada fold, o pipeline é ajustado **apenas nos dados de treino**; scaler e PCA são recalculados a cada iteração, evitando vazamento de informação.
  - Métricas por fold (acurácia, sensibilidade, especificidade, AUC, matriz de confusão) são registradas em `TABELAS/resultados_cv_regressao_logistica.csv`.
  - Em paralelo, usa‑se `cross_val_predict(..., method="predict_proba")` para obter probabilidades out‑of‑fold e calcular uma **AUC global out‑of‑fold** e uma matriz de confusão global. As figuras correspondentes são:
    - `FIGURAS/matriz_confusao_regressao_logistica.png`;
    - `FIGURAS/roc_regressao_logistica.png`.

- **Interpretação**
  - O pipeline é ajustado em todo o conjunto `X_log/Y` para fins descritivos, e os coeficientes da regressão logística no espaço dos PCs são salvos em `TABELAS/coeficientes_regressao_logistica_pcs.csv`.
  - Uma importância aproximada por gene é derivada combinando esses coeficientes com os loadings do PCA do pipeline (`TABELAS/importancia_genes_modelo.csv`), explicitamente como análise auxiliar, não como coeficientes genuínos no espaço gênico original.

### Resultados e discussão

- **PCA exploratório**  
  O gráfico PC1 vs PC2 (`pca_mut_vs_gc.png`) mostra agrupamentos consistentes por tipo celular e condição, com **separação parcial** entre amostras mut e GC. Há clusters em que mut e GC se separam de forma razoável e outros em que permanecem misturados, refletindo a forte influência da estrutura experimental (iPSC/MSC/NSC; Ctrl/UV/RS).

- **Desempenho do modelo**  
  A regressão logística com PCA (5 PCs) apresentou:
  - **acurácia média em CV** ≈ 0,76;
  - **sensibilidade média** (mut como classe positiva) ≈ 0,63;
  - **especificidade média** ≈ 0,87;
  - **AUC média por fold** ≈ 0,73.
  A **AUC global out‑of‑fold**, calculada a partir das probabilidades produzidas pelo pipeline em `cross_val_predict`, foi ≈ 0,74, consistente com a curva ROC (`roc_regressao_logistica.png`).

- **Matriz de confusão global**  
  A matriz de confusão out‑of‑fold (`matriz_confusao_regressao_logistica.png`) indica:
  - 9 controles GC corretamente classificados e 3 rotulados como mut;
  - 8 mut corretamente classificados e 4 rotulados como GC.  
  Isso corresponde a uma acurácia global em torno de 0,7, com **boa especificidade** e **sensibilidade moderada**.

- **Interpretação dos coeficientes**  
  Os coeficientes da regressão logística nos PCs mostram maior peso em um dos componentes (PC5), sugerindo que parte da separação mut vs GC está concentrada em uma direção específica no espaço reduzido. A análise de importância aproximada por gene fornece um ranking útil para exploração biológica, mas foi tratada explicitamente como **heurística** e não como estimativa de efeito causal.

- **Limitações**  
  - Tamanho amostral reduzido (24 amostras), o que torna as estimativas de desempenho instáveis; alguns folds apresentam AUC próxima de 1, outros bem menores.  
  - Estrutura experimental complexa, com diferentes tipos celulares e condições, e combinações incompletas entre eles. Isso dificulta o uso de esquemas de validação por grupos biológicos mais rigorosos e reforça o caráter **exploratório** do modelo.  
  - Uso de contagens NCBI‑generated (e não necessariamente idênticas às do artigo original), o que implica que comparações diretas com resultados publicados devem ser feitas com cautela.

### Conclusão

O pipeline implementado mostra que há um **sinal transcriptômico detectável** associado ao estado mutante vs gene‑corrected em células derivadas de paciente com síndrome de Cockayne no conjunto GSE124208. A regressão logística em um espaço de 5 componentes principais alcança desempenho moderado (AUC global ≈ 0,74), com maior segurança em identificar controles GC do que mutantes.  
Esses resultados apoiam a interpretação de que o estado mut/GC se reflete em padrões globais de expressão gênica, mas, dadas as limitações de tamanho amostral e de desenho experimental, o modelo deve ser entendido como **ferramenta exploratória de separação de classes**, não como classificador diagnóstico robusto. O código e as saídas foram organizados para facilitar reanálises futuras, ajustes de hiperparâmetros (por exemplo, variação no número de PCs) e extensões específicas por tipo celular.