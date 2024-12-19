# Анализ данных 16S

## Оценка качества и фильтрация прочтений

Нам было предоставлено 36 файлов в формате FASTA, а также метаданные, включающие информацию о группе образцов (sample 1 или 2), использовании антибиотиков и дне взятия пробы.

### Оценка качества

Качество прочтений было оценено с использованием FastQC и суммировано с помощью MultiQC. На основании анализа было установлено, что качество прочтений падает с 85 по 130 нуклеотид. Это потребовало проведения фильтрации данных.

![График качества прочтений](pic/pic_fastqc.png)

### Фильтрация по качеству

Для фильтрации данных использовался инструмент Trimmomatic с параметрами:

- **SLIDINGWINDOW:4:20** — обрезка по качеству с использованием скользящего окна длиной 4 и порогом качества 20.
- **MINLEN:50** — удаление прочтений длиной менее 50 нуклеотидов.

Настройки фильтрации были реализованы в скрипте cut.sh.

Дальнейший анализ проводился в R с использованием пакета DADA2.

Данные были дополнительно отфильтрованы с помощью функции `filterAndTrim`. Были использованы следующие параметры:

- **maxN=0** — удаляются прочтения с любым количеством неидентифицированных нуклеотидов (N).
- **maxEE=2** — допускается не более двух ожидаемых ошибок для прочтения.
- **truncQ=2** — обрезка прочтений при качестве ниже 2.
- **rm.phix=TRUE** — удаление последовательностей PhiX (контрольная последовательность для Illumina).

### Ординация данных

Для визуализации структуры данных была проведена ординация с использованием метода главных компонент (PCA) на основе Евклидовой метрики. 

![PCA на основе Евклидовой метрики](pic/PCA_euclid.png)

Результаты PCA с использованием двух компонент не выявили чёткого разделения образцов на кластеры. В связи с этим было решено использовать другие метрики схожести: **косинусное расстояние** и **сходство Брея-Кертиса**.

Эти методы показали визуальное разделение образцов на три кластера:

<div style="display: flex; justify-content: space-between;">
    <img src="pic/PCA_cos.png" alt="PCA на основе косинусного расстояния" width="45%">
    <img src="pic/PCA_bray.png" alt="PCA на основе сходства Брея-Кертиса" width="45%">
</div>

Для дальнейшего анализа было решено использовать **сходство Брея-Кертиса**.

PCA-плот с отображением метаданных показал, что разделение на кластеры **не зависит** от следующих факторов:
- Донор,
- Приём антибиотиков,
- Дата взятия пробы.

![PCA с отображением метаданных](pic/PCA_w_metadata.png)

### Кластеризация

Поскольку разделение на кластеры не связано с указанными метаданными, было предположено наличие других факторов. В данные был добавлен номер кластера, и кластеризация была проведена с использованием метода **k-means**.

### Анализ таксономии

Для каждого кластера были выделены 5 наиболее часто встречающихся таксонов. Полный список представлен ниже:

| Кластер | Таксономия                                                                                                      | Общая встречаемость |
|---------|----------------------------------------------------------------------------------------------------------------|---------------------:|
| 1       | Bacteria; Bacillota; Bacilli; Staphylococcales; Staphylococcaceae; Staphylococcus; aureus                      |               17779 |
| 1       | Bacteria; Bacillota; Bacilli; Lactobacillales; Streptococcaceae; Streptococcus                                 |               11541 |
| 1       | Bacteria; Pseudomonadota; Gammaproteobacteria; Pseudomonadales; Pseudomonadaceae; Pseudomonas; tolaasii        |                5409 |
| 1       | Bacteria; Cyanobacteriota; Cyanobacteriia; Chloroplast                                                         |                4930 |
| 1       | Bacteria; Actinomycetota; Actinobacteria; Mycobacteriales; Corynebacteriaceae; Corynebacterium; tuberculostearicum |             4004 |
| 2       | Bacteria; Bacteroidota; Bacteroidia; Bacteroidales; Bacteroidaceae; Bacteroides; vulgatus                      |               12488 |
| 2       | Bacteria; Bacteroidota; Bacteroidia; Bacteroidales; Bacteroidaceae; Bacteroides; dorei                         |                9138 |
| 2       | Bacteria; Bacteroidota; Bacteroidia; Bacteroidales; Bacteroidaceae; Bacteroides; ovatus                        |                9083 |
| 2       | Bacteria; Verrucomicrobiota; Verrucomicrobiia; Verrucomicrobiales; Akkermansiaceae; Akkermansia; muciniphila   |                6839 |
| 2       | Bacteria; Bacteroidota; Bacteroidia; Bacteroidales; Bacteroidaceae; Bacteroides; uniformis                     |                6826 |
| 3       | Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Pasteurellaceae; Haemophilus                  |                9987 |
| 3       | Bacteria; Pseudomonadota; Gammaproteobacteria; Burkholderiales; Neisseriaceae; Neisseria; perflava             |                8826 |
| 3       | Bacteria; Bacillota; Bacilli; Lactobacillales; Streptococcaceae; Streptococcus                                 |                5144 |
| 3       | Bacteria; Fusobacteriota; Fusobacteriia; Fusobacteriales; Fusobacteriaceae; Fusobacterium; periodonticum       |                4744 |
| 3       | Bacteria; Bacteroidota; Bacteroidia; Bacteroidales; Prevotellaceae; Prevotella; melaninogenica                 |                4493 |

На основании анализа таксономии было определено соответствие кластеров предполагаемой локации взятия пробы:

- **Кластер 1**: Бактерии, характерные для кожи человека.
- **Кластер 2**: Бактерии, характерные для кишечника.
- **Кластер 3**: Бактерии, характерные для ротовой полости.

Эти результаты позволяют предполагать, что разделение на кластеры связано с происхождением проб из разных частей тела человека.
