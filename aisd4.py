#20. Формируется матрица F следующим образом: скопировать в нее А и  если в Е количество чисел, больших К в четных столбцах ,
# чем произведение чисел в нечетных строках , то поменять местами С и Е симметрично, иначе С и В поменять местами несимметрично. При этом матрица А не меняется.
# 1После чего если определитель матрицы А больше суммы диагональных элементов матрицы F, то вычисляется выражение: A-1*AT – K * FТ,
# иначе вычисляется выражение (A +GТ-F-1)*K, где G-нижняя треугольная матрица, полученная из А. Выводятся по мере формирования А, F и все матричные операции последовательно.

import matplotlib.pyplot as mpl #Импорт библиотек
import numpy as np

K_test = 3 # Тестовые данные
N_test = 10
A_test = np.array([[9, 5, 1, 6, -3, -8, -7, 1, 1, 10],
    [1, 6, -5, -1, -4, -1, 10, 5, -10, -6],
    [-8, -2, -3, 7, 9, 1, 8, 0, 9, 5],
    [-7, 6, 0, -8, 4, 2, 1, -8, -5, -1],
    [-3, -9, -4, -1, -5, -3, -6, 9, 7, -6],
    [-7, 1, 7, 8, -3, 5, 7, -1, -7, -6],
    [-1, 6, -5, 2, 2, 2, 3, 10, -8, 4],
    [-4, -2, 1, -2, -2, -4, -7, -10, 15, 5],
    [2, -3, 0, -7, -1, 0, 9, -8, 9, 4],
    [-8, -10, 3, 0, -5, 10, -8, -10, -1, 8]
    ])

print('Использовать тестовые данные или случайные?')
while True:
    choice = input('Введите 1, если хотите использовать тестовые данные, 2 - если случайные, q - для выхода из программы): ')
    if choice == '1' or choice == '2' or choice == 'q':
        break

if choice == '1':
    K = K_test
    N = N_test
    A = A_test

if choice == '2': # Генерация случайных данных
    K = int(input("Введите число К="))
    while True:
        N = int(input("Введите число N="))
        if N < 6:
            print('Число N слишком малое. Введите N >= 6')
        else:
            break
# Формируем матрицу А
    A = np.random.randint(-10, 10, size=(N, N))

if choice == 'q':
    exit()

n = N // 2  # Размерность матриц B, C, D, E (n x n)
n_first = n
if N % 2 == 0:
    n_second = n
else:
    n_second = n+1

E = A[:n_first,:n_first]
B = A[:n_first,n_second::]
C = A[n_second::,n_second::]
D = A[n_second::,:n_first]

# Печатаем матрицы A, E, B, C, D
print('\nМатрица A:\n', A)
print('\nМатрица E:\n', E)
print('\nМатрица B:\n', B)
print('\nМатрица C:\n', C)
print('\nМатрица D:\n', D)

# Формируем матрицу F
count = 0
mult = 1
for i in E[::,1::2]:
    for j in i:
        if j > K:
            count +=1

for i in E[::2]:
    for j in i:
        mult *= j

if count < mult:
    F = np.hstack((A[::,:n_second], A[::-1,n_second::]))
    print('\nМатрица F, подматрицы C и E симметрично заменены:\n', F)
else:
    x = A[::,n_second::]
    for i in x:
        for j in range(n_second):
            x[j], x[N-1-j] = x[N-1-j], x[j].copy()
    F = np.hstack((A[::,:n_second], x))
    print('\nМатрица F, подматрицы C и B несимметрично заменены:\n', F)

# Вычисляем выражение в зависимости от условия
if np.linalg.det(A) > np.trace(F):
    print('\nВычисляем выражение A-1*AT – K * FТ')
    inv_A = np.linalg.inv(A)
    print('\nОбратная матрица A:\n', inv_A)
    A_trans = A.transpose()
    print('\nТранспонированная матрица A:\n', A_trans)
    inv_A_and_A_trans = np.dot(inv_A, A_trans)
    print('\nУмножение обратной и транспонированной матрицы A:\n',1, inv_A_and_A_trans)
    F_trans = F.transpose()
    print('\nТранспонированная матрица F:\n', F_trans)
    K_and_F_trans = np.dot(K, F_trans)
    print('\nУмножение константы К и транспонированной матрицы F:\n', K_and_F_trans)
    Result = np.subtract(inv_A_and_A_trans, K_and_F_trans)
else:
    print('\nВычисляем выражение (A +GТ-F-1)*K')
    A_trans = A.transpose()
    print('\nТранспонированная матрица A:\n', A_trans)
    G = np.tril(A)
    print('\nНижняя треугольная матрица G из матрицы A:\n', G)
    F_trans = F.transpose()
    print('\nТранспонированная матрица F:\n', F_trans)
    A_trans_plus_G = np.add(A_trans, G)
    print('\nСумма транспонированной матрицы А и треугольной матрицы G:\n', G)
    A_t_plus_G_minus_F_t = np.subtract(A_trans_plus_G, F_trans)
    print('\nРазность предыдущего выражения и транспонированной матрицы F:\n', G)
    Result = np.dot(A_t_plus_G_minus_F_t, K)
print('\nРезультирующая матрица:\n', Result)

mpl.figure(figsize=(16, 9))

# тепловая карта матрицы F
mpl.subplot(1, 3, 1)
mpl.xticks(ticks=np.arange(F.shape[1]))
mpl.yticks(ticks=np.arange(F.shape[1]))
mpl.xlabel('Номер столбца')
mpl.ylabel('Номер строки')
hm = mpl.imshow(F, cmap='plasma', interpolation="nearest")
mpl.colorbar(hm, shrink = 0.5)
mpl.title('Тепловая карта элементов')

#График максимальных элементов в столбцах матрицы F
res = np.amax(F, axis=0)
x = np.arange(F.shape[1])
mpl.subplot(1, 3, 2)
mpl.bar(x, res, label='Максимальный элемент в столбце')
mpl.xlabel('Номер столбца')
mpl.ylabel('Максимальный элемент')
mpl.title('График максимальных элементов в столбцах')
mpl.legend()

# круговая диаграмма матрицы F
x = np.arange(F.shape[1])
mpl.subplot(1, 3, 3)
P = []
for i in range(N):
    P.append(abs(F[0][i]))
mpl.pie(P, labels=x, autopct='%1.2f%%')
mpl.title("Диаграмма pie")

mpl.tight_layout(pad=2.5, w_pad=1, h_pad=1) # расстояние от границ и между областями
mpl.suptitle("Графический вывод", y=1)
mpl.show()
