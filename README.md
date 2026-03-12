Краткая инфа

Что реализовано
1) Разбор FASTQ
- Функция: iter_fastq(path)
- читает файл по 4 строкам FASTQ и возвращает пары `(sequence, quality)`

2) Функция: task1(path)
- считает число ридов, минимальную, среднюю (с округлением) и максимальную длину

3) Функция: task2(path)
- считает процент G и C по всем нуклеотидам

4) Функция: task3(path, position=10)
- считает средний `Phred` в фиксированной позиции (1-based)

5) Тримминг скользящим окном и фильтр по длине
- Функция: `slidingwindow_len(sequence, quality, window=5, required_quality=30)`
- тримминг одного рида по логике `SLIDINGWINDOW` (как в Trimmomatic)

- Функция: `task4(path, window=5, quality=30, min_len=60)`
- применяет тримминг ко всем ридам и возвращает:
  - число удаленных ридов после первого шага,
  - минимальную/среднюю/максимальную длину в отфильтрованном файле,
  - число ридов после фильтра `MINLEN >= 60`

Чем проверяли
- Тесты
- Сверка с Trimmomatic:
  - java -jar trimmomatic-0.39.jar SE -phred33 reads.fastq.txt out.fastq SLIDINGWINDOW:5:30 MINLEN:60
