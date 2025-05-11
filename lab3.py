#Цели
##Освоить методы визуализации данных с использованием диаграммы рассеяния. 
##Научиться разделять данные на классы с помощью цветового кодирования. Практиковать работу с библиотеками для анализа и визуализации данных.
#Инструменты
##Библиотеки matplotlib и seaborn для визуализации. Набор данных iris загружен через scikit-learn.
#Выводы
##Диаграмма рассеяния отображает распределение данных по двум признакам. 
##Цветовое кодирование позволяет четко выделить три класса ирисов. 
##Инструменты визуализации (seaborn, matplotlib) эффективны для анализа многоклассовых данных.

import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
import seaborn as sns


iris = load_iris()
data = iris.data
target = iris.target
feature_names = iris.feature_names


x = data[:, 0]  # sepal length
y = data[:, 1]  # sepal width


plt.figure(figsize=(10, 6))
sns.scatterplot(x=x, y=y, hue=target, palette='viridis', s=70, edgecolor='k')
plt.title('Диаграмма рассеяния для набора Iris')
plt.xlabel(feature_names[0])
plt.ylabel(feature_names[1])
plt.legend(title='Классы', labels=iris.target_names)
plt.grid(True)
plt.show()