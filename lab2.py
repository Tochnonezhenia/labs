#Цели  
##Освоить методы анализа GC-состава последовательностей ДНК. 
##Реализовать алгоритм обработки FASTA-файлов с использованием библиотеки Biopython. 
##Обеспечить надежность работы программы через проверку существования файла, прав доступа и обработку ошибок. 
##Определить последовательность с максимальным значением GC-состава среди всех записей.  

#Инструменты  
##Библиотека Biopython (модуль SeqIO) для чтения и парсинга FASTA-файлов.
##Системные функции (os.path.exists, os.access) для проверки доступности файла и прав доступа. 
##Механизмы обработки ошибок (try-except, sys.stderr.write) для повышения устойчивости программы.  

#Выводы  
##Программа корректно вычисляет GC-состав для последовательностей ДНК из FASTA-файла. 
##Реализованы проверки на существование файла, его формат и права доступа, что исключает частые ошибки ввода. 
##Алгоритм эффективно обрабатывает данные длиной до 1 килобазы с линейной сложностью. 
##Ограничения включают отсутствие валидации символов (A, T, C, G) и округления результатов. 
##Для улучшения рекомендуется добавить проверку на некорректные символы и расширить функционал для анализа нескольких файлов.
from Bio import SeqIO
import os
import sys
#указываем свой файл
fastafile = "seq.fasta"

#существования файла
if not os.path.exists(fastafile):
    sys.stderr.write(f"Ошибка: Файл '{fastafile}' не найден.\n")
    sys.exit(1)

#это именно файл
if not os.path.isfile(fastafile):
    sys.stderr.write(f"Ошибка: '{fastafile}' не является файлом.\n")
    sys.exit(1)

#права на чтение
if not os.access(fastafile, os.R_OK):
    sys.stderr.write(f"Ошибка: Нет прав на чтение файла '{fastafile}'.\n")
    sys.exit(1)

#права на удаление (права на запись в директорию)
dir_name = os.path.dirname(fastafile) or '.'  # Текущая директория если путь отсутствует
if not os.access(dir_name, os.W_OK):
    sys.stderr.write(f"Ошибка: Нет прав на удаление файла в директории '{dir_name}'.\n")
    sys.exit(1)

try:
    records = list(SeqIO.parse(fastafile, "fasta"))
except Exception as e:
    sys.stderr.write(f"Ошибка при чтении файла: {str(e)}\n")
    sys.exit(1)

res = {}
for x in records:
    GC_count = (x.count('C') + x.count('G')) / len(x) * 100
    res[x.id] = GC_count

if res:
    max_procent = max(res, key=res.get)
    print(max_procent, res[max_procent])
else:
    sys.stderr.write("Файл не содержит корректных FASTA записей.\n")
    sys.exit(1)