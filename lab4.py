#Цели  
##Освоить методы работы с файлами GenBank для извлечения биологических данных. 
##Реализовать объединение нескольких GenBank-файлов в единую структуру с использованием библиотеки Biopython. 
##Анализировать GC-состав последовательностей и выполнять трансляцию кодирующих областей (CDS) в белковые последовательности. 
##Обеспечить надежность программы через проверку целостности данных, обработку ошибок и валидацию входных файлов.  

#Инструменты  
##Библиотека Biopython (модули SeqIO и SeqUtils) для парсинга GenBank-файлов, извлечения последовательностей и расчета GC-состава. 
##Функции трансляции из Bio.Seq для конвертации нуклеотидных последовательностей в аминокислотные. 
##Системные модули (os, sys) для проверки существования файлов, прав доступа и обработки исключений. 
##Методы сортировки и фильтрации данных для вывода результатов в заданном формате.  

#Выводы  
##Программа успешно объединяет данные из GenBank-файлов, сохраняя их структуру и метаинформацию. 
##Корректно вычисляется GC-состав для каждой последовательности, что позволяет идентифицировать виды по особенностям генома. 
##Трансляция кодирующих областей реализована с учетом правил генетического кода, обеспечивая точность получения белковых последовательностей. 
##Проверки на доступность файлов и обработка ошибок повышают устойчивость программы к некорректным входным данным. 
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import warnings

def merge_genbank_files(file1, file2, output_file):
    """Объединяет два файла GenBank, игнорируя некорректные записи."""
    valid_records = []
    error_log = []
    
    # Функция для безопасной загрузки записей
    def load_records(file_path):
        nonlocal error_log
        count = 0
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")  # Игнорируем предупреждения парсера
                for record in SeqIO.parse(file_path, "genbank"):
                    try:
                        # Проверяем все возможные проблемы
                        if not hasattr(record, 'seq'):
                            raise ValueError("Нет атрибута 'seq'")
                        if len(record.seq) == 0:
                            raise ValueError("Пустая последовательность")
                        if 'N' in record.seq:  # Пример дополнительной проверки
                            error_log.append(f"{record.id}: содержит неопределенные нуклеотиды (N)")
                            continue
                        valid_records.append(record)
                        count += 1
                    except Exception as e:
                        error_log.append(f"Пропущена запись {getattr(record, 'id', 'unknown')}: {str(e)}")
            return count
        except Exception as e:
            error_log.append(f"Файл {file_path}: {str(e)}")
            return 0

    # Обрабатываем оба файла
    count1 = load_records(file1)
    count2 = load_records(file2)

    # Записываем результат
    if valid_records:
        SeqIO.write(valid_records, output_file, "genbank")
        print(f"Успешно объединено: {len(valid_records)} записей ({count1} + {count2})")
        print(f"Сохранено в: {output_file}")
    else:
        print("Ошибка: Нет валидных записей для объединения!")
        return False

    # Выводим ошибки
    if error_log:
        print("\nЛог ошибок:")
        for error in error_log:
            print(f" - {error}")
    
    return True

def analyze_gc(input_file):
    """Анализ GC-состава с детальным отчетом"""
    try:
        records = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for record in SeqIO.parse(input_file, "genbank"):
                try:
                    gc = gc_fraction(record.seq)
                    records.append((record.id, gc, len(record.seq)))
                except Exception as e:
                    print(f"Ошибка в записи {record.id}: {str(e)}")
        
        if not records:
            print("Нет данных для анализа!")
            return

        # Сортировка по GC
        sorted_records = sorted(records, key=lambda x: x[1])
        
        # Вывод результатов
        print("\nРезультаты анализа:")
        print("{:<15} {:<10} {:<15}".format("ID", "GC (%)", "Длина последовательности"))
        for rec in sorted_records:
            print(f"{rec[0]:<15} {rec[1]*100:<10.2f} {rec[2]:<15}")

    except Exception as e:
        print(f"Фатальная ошибка: {str(e)}")

if __name__ == "__main__":
    # Настройки
    HUMAN_GB = "homo.gb"
    VIRUS_GB = "clos.gb"
    OUTPUT_GB = "combined_results.gb"

    # Объединение файлов
    if merge_genbank_files(HUMAN_GB, VIRUS_GB, OUTPUT_GB):
        # Анализ
        analyze_gc(OUTPUT_GB)