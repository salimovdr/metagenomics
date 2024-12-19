# Метагеномный анализ

Для анализа был выбран образец TB2_SN3. Основной целью было определить таксономический состав образца и выполнить биннинг для получения высококачественных MAG.

## Подготовка данных
- **Исходный файл:** `scaffolds_TB2_SN3.fasta`
- Данные были подготовлены для анализа с использованием инструментов `kraken2` для таксономической классификации и `metabat2`/`maxbin2` для биннинга.

**Результаты:**
- Классифицированные последовательности: 12,20%
- Неклассифицированные последовательности: 87,81%

Отчет был визуализирован с помощью инструмента `Pavian`.

---

### 3. Биннинг
Процесс биннинга был выполнен с использованием двух инструментов:

#### MetaBAT2
**Параметры запуска:**
- Входной файл: `scaffolds_TB2_SN3.fasta`
- Выходная папка: `metabat_results/`
- Потоки: 10

#### MaxBin2
**Параметры запуска:**
- Входной файл: `scaffolds_TB2_SN3.fasta`
- Покрытие: рассчитано из названий контигов
- Выходная папка: `maxbin_results/`
- Потоки: 10

**Результаты:**
- MetaBAT2 сгенерировал несколько бинов с разной полнотой и контаминацией.
- MaxBin2 также сгенерировал бины, но ни один из них не соответствовал критериям MIMAG.

---

### 4. Оценка качества бинов
Качество бинов оценивалось с использованием инструмента `CheckM`.

**Параметры запуска:**
- Входная папка: `metabat_results/`
- Выходная папка: `checkm_metabat_results/`
- Расширение файлов: `fa`
- Потоки: 10

#### Результаты:
- Один бин (`bin.2`) из MetaBAT2 соответствует критериям MIMAG для **чернового генома среднего качества**:
  - Полнота: 58,88%
  - Контаминация: 8,81%
  
Бины из MaxBin2 не соответствуют необходимым стандартам качества.

---

## Оставшиеся шаги
1. **Аннотация:** Использовать `Prokka` для аннотации выбранного бина.
2. **Анализ кластеров:** Выполнить анализ с использованием AntiSMASH для поиска генов вторичных метаболитов.
3. **Завершение отчета:** Дополнить результаты после выполнения оставшихся этапов.

---

## Заключение
В ходе анализа был успешно идентифицирован один MAG среднего качества из образца `TB2_SN3`. Следующие шаги включают аннотацию и изучение потенциальных кластеров вторичных метаболитов.