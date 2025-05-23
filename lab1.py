#Цели
##Освоить метод подсчёта расстояния Хэмминга между двумя строками ДНК. 
##Реализовать алгоритм для определения минимального числа точечных мутаций. 
##Обеспечить проверку корректности входных данных (равенство длин строк).
#Инструменты
##Функция count_mismatches для посимвольного сравнения строк. Обработка ввода через input().
#Выводы
##Разработанная программа корректно вычисляет расстояние Хэмминга между двумя строками ДНК. 
##Реализована проверка на равенство длин входных данных. Алгоритм демонстрирует эффективность для строк длиной до 1000 символов. 
##Метод может быть применён для анализа эволюционных изменений в биологических последовательностях.
def count_mismatches(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Строки должны быть одинаковой длины")
    mismatches = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            mismatches += 1
    return mismatches

string1 = input("Введите первую строку: ")
string2 = input("Введите вторую строку: ")
print(f"Количество несовпадающих элементов: {count_mismatches(string1, string2)}")
