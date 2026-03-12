from pathlib import Path

#Итерирует FASTQ-файл и возвращает пары (sequence, quality) по одному прочтению
def iter_fastq(path: Path):
    with path.open("r", encoding="utf-8") as f:
        while True:
            if not f.readline():  
                break
            seq = f.readline().strip()
            f.readline()  
            quality = f.readline().strip()
            yield seq, quality

#Задание 1: базовая статистика длин прочтений
def task1(path: Path) -> tuple[int, int, int, int]:
    lengths = [len(seq) for seq, _ in iter_fastq(path)]
    total_reads = len(lengths)
    min_len = min(lengths)
    avg_len = round(sum(lengths) / total_reads)
    max_len = max(lengths)
    return total_reads, min_len, avg_len, max_len

#Задание 2: GC-состав по всем нуклеотидам файла в процентах
def task2(path: Path) -> float:
    gc = 0
    total = 0
    for seq, _ in iter_fastq(path):
        upper = seq.upper()
        gc += upper.count("G") + upper.count("C")
        total += len(upper)
    return round(gc * 100 / total, 2)

#Задание 3: среднее качество Phred в заданной позиции, округленное до целого
def task3(path: Path, position: int = 10) -> int:
    qsum = 0
    count = 0
    idx = position - 1
    for seq, qual in iter_fastq(path):
        if idx < len(seq):
            qsum += ord(qual[idx]) - 33
            count += 1
    return round(qsum / count)


#Задание 4: применяет SLIDINGWINDOW ко всем ридам и считает статистику до/после MINLEN
def task4(path: Path, window: int = 5, quality: int = 30, min_len: int = 60) -> dict:
    total_reads = 0
    dropped_reads = 0
    lengths_after_sw = []
    for sequence, qline in iter_fastq(path):
        total_reads += 1
        trimmed_len = slidingwindow_len(sequence, qline, window=window, required_quality=quality)
        if trimmed_len is None:
            dropped_reads += 1
            continue
        lengths_after_sw.append(trimmed_len)
    min_length = min(lengths_after_sw)
    mid_length = round(sum(lengths_after_sw) / len(lengths_after_sw))
    max_length = max(lengths_after_sw)
    reads_after_minlen = sum(1 for length in lengths_after_sw if length >= min_len)
    
    return {
        "dropped_reads": dropped_reads,
        "min_len": min_length,
        "mid_len": mid_length,
        "max_len": max_length,
        "after_minlen": reads_after_minlen,
    }


#Тримминг одного прочтения; возвращает новую длину 
def slidingwindow_len(sequence: str, quality: str, window: int = 5, required_quality: int = 30) -> int | None:
    quals = [
        0 if base == "N" else (ord(qchar) - 33)
        for base, qchar in zip(sequence, quality)
    ]
    if len(quals) < window:
        return None
    required_total = required_quality * window
    total = sum(quals[:window])
    if total < required_total:
        return None
    len_to_keep = len(quals)
    for i in range(0, len(quals) - window):
        total = total - quals[i] + quals[i + window]
        if total < required_total:
            len_to_keep = i + window
            break
    i = len_to_keep
    while i > 1 and quals[i - 1] < required_quality:
        i -= 1
    if i < 1:
        return None
    return i


path = Path("reads.fastq.txt")
total, min_len, mid_len, max_len = task1(path)
gc = task2(path)
q10 = task3(path, position=10)
t4 = task4(path, window=5, quality=30, min_len=60)

print("Task 1")
print(f"Reads: {total}")
print(f"Min len: {min_len}")
print(f"mid len: {mid_len}")
print(f"Max len: {max_len}")
print()
print("Task 2")
print(f"GC: {gc}")
print()
print("Task 3")
print(f"Mean phred at position 10: {q10}")
print()
print("Task 4")
print(f"Dropped only: {t4['dropped_reads']}")
print(f"Min/Avg/Max after SW: {t4['min_len']} / {t4['mid_len']} / {t4['max_len']}")
print(f"After MINLEN:60: {t4['after_minlen']}")


