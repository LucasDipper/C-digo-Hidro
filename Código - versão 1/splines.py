import matplotlib.pyplot as plt   # biblioteca de gráficos (importada, ainda não usada no código)
import numpy as np                # usado para np.linspace, np.array e np.cross

# Tabela de cotas de EXEMPLO (inventada por nós, não é dado real de embarcação).
# Formato de cada item: (x_baliza, [alturas z conhecidas], [meias-bocas y conhecidas])
# Repare que popa e proa têm alcance menor que a meia-nau - isso imita a irregularidade
# de uma tabela de cotas real (perto das pontas, o casco não alcança as linhas d'água mais altas).
malha_teste = [
    (0.0,  [0.0, 2.0],       [0.0, 1.0]),        # popa - só tem dado até z=2
    (5.0,  [0.0, 2.0, 4.0],  [0.0, 2.0, 3.0]),    # meia-nau - tem dado até z=4 (o convés)
    (10.0, [0.0, 2.0],       [0.0, 1.0]),         # proa - simétrica à popa
]

# ============================================================
# MOTOR DE SPLINE CÚBICA NATURAL (curva única: um x, um y)
# ============================================================

def calcular_h_num(x):
    # Calcula os espaçamentos entre pontos consecutivos de UMA curva,
    # e devolve também n e num_incognitas já prontos, para não recalcular
    # esses dois valores em cada uma das próximas funções.
    n = len(x)                       # quantidade de pontos conhecidos da curva
    h = []                           # vai guardar os n-1 espaçamentos entre pontos
    num_incognitas = n - 2           # nº de g's a resolver (spline natural: g0 = g(n-1) = 0, ficam de fora)
    for i in range(n - 1):           # n pontos -> n-1 intervalos
        h.append(x[i + 1] - x[i])    # espaçamento entre o ponto i e o ponto i+1
    return h, n, num_incognitas      # devolve os três de uma vez, prontos para as próximas funções

def montar_sistema(y, h, num_incognitas):
    # Monta a matriz tridiagonal e o vetor do sistema linear A*g = vetor,
    # cuja solução são os g_k (segundas derivadas da spline em cada nó).
    matriz = [[0] * num_incognitas for _ in range(num_incognitas)]  # matriz quadrada, começa zerada
    vetor  = [0] * num_incognitas                                   # lado direito do sistema, começa zerado

    for linha in range(num_incognitas):     # uma linha do sistema para cada incógnita g_k
        k = linha + 1                        # k percorre 1, ..., n-2 (índice real do g, já que g0/g(n-1) não entram)

        # coeficientes da equação geral:
        # h[k-1]*g[k-1] + 2*(h[k-1]+h[k])*g[k] + h[k]*g[k+1] = lado_direito
        coef_esq  = h[k - 1]                        # multiplica g(k-1)
        coef_meio = 2 * (h[k - 1] + h[k])            # multiplica g(k) - sempre existe, é a diagonal principal
        coef_dir  = h[k]                             # multiplica g(k+1)

        # posiciona os coeficientes na linha certa da matriz
        if linha - 1 >= 0:
            matriz[linha][linha - 1] = coef_esq      # só existe se não for a primeira linha
        matriz[linha][linha] = coef_meio             # diagonal principal, sempre presente
        if linha + 1 < num_incognitas:
            matriz[linha][linha + 1] = coef_dir      # só existe se não for a última linha

        # lado direito: 6 * ( (y[k+1]-y[k])/h[k] - (y[k]-y[k-1])/h[k-1] )
        vetor[linha] = 6 * ((y[k + 1] - y[k]) / h[k] - (y[k] - y[k - 1]) / h[k - 1])

    return matriz, vetor   # sistema pronto para ser resolvido pela eliminação de Gauss

def eliminacao_gauss(matriz, vetor, num_incognitas):
    # Resolve o sistema linear (matriz, vetor) por eliminação de Gauss:
    # primeiro escalona (zera abaixo da diagonal), depois resolve de baixo para cima.

    # --- ESCALONAMENTO ---
    for coluna in range(num_incognitas - 1):                # percorre cada coluna a ser zerada
        for l_alvo in range(coluna + 1, num_incognitas):     # percorre as linhas abaixo da diagonal nessa coluna
            fator = matriz[l_alvo][coluna] / matriz[coluna][coluna]  # fator que zera o elemento
            vetor[l_alvo] -= fator * vetor[coluna]            # aplica a mesma operação no vetor

            for i in range(coluna, num_incognitas):           # aplica a subtração em toda a linha da matriz
                matriz[l_alvo][i] -= fator * matriz[coluna][i]

    # --- SUBSTITUIÇÃO REVERSA ---
    g = [None] * num_incognitas                                       # vai guardar a solução (os g_k)
    for i in range(num_incognitas - 1, -1, -1):  # começa pelo mesmo numero de incógnitas, vai até 0, diminuindo 1
        soma = 0
        for j in range(num_incognitas - 1, -1, -1):
            if j > i:                              # soma as incógnitas à direita de i, já resolvidas antes
                soma += g[j]*matriz[i][j]
        g[i] = (vetor[i] - soma)/matriz[i][i]       # isola g[i] na equação

    return matriz, vetor, g   # matriz/vetor já escalonados (não usados depois) + a solução g

def calcula_coeficientes(y, h, g, n):
    # Calcula os coeficientes a, b, c, d de cada trecho da spline (Eq. 3.3
    # do relatório), a partir dos g_k já resolvidos.
    g_completo = [0] + g + [0]   # reinsere g0=0 e g(n-1)=0, removidos do sistema pela condição natural
    a = [None] * n               # coeficiente do termo cúbico de cada trecho
    b = [None] * n               # coeficiente do termo quadrático
    c = [None] * n               # coeficiente do termo linear
    d = [None] * n               # termo constante (valor no nó da direita do trecho)

    for k in range(1, n):        # um conjunto de coeficientes para cada trecho (k = 1 até n-1)
        a[k] = (g_completo[k] - g_completo[k-1])/(6*h[k-1])
        b[k] = g_completo[k]/2
        c[k] = ((y[k] - y[k-1])/h[k-1]) + ((2*h[k-1]*g_completo[k] + g_completo[k-1]*h[k-1])/6)
        d[k] = y[k]               # o termo constante é sempre o valor de y no nó da direita

    return a, b, c, d

def avalia_spline(x_alvo, x, y, a, b, c, d):
    # Avalia a spline já ajustada num ponto x_alvo qualquer: primeiro acha
    # em qual trecho [x[k-1], x[k]] ele cai, depois calcula o polinômio ali.
    n = len(x)
    for k in range(1, n):                 # percorre os nós procurando o trecho certo
        if x[k-1] <= x_alvo <= x[k]:
            break                          # k fica com o índice do trecho encontrado

    x_k = x[k]   # nó da direita do trecho (referência usada na fórmula do polinômio)
    resultado = a[k]*(x_alvo - x_k)**3 + b[k]*(x_alvo - x_k)**2 + c[k]*(x_alvo - x_k) + d[k]
    return resultado

def ajustar_spline(x, y):
    # Encapsula a sequência inteira: resolve o sistema UMA VEZ para a curva
    # (x, y) e devolve os coeficientes prontos, para avaliar quantas vezes quiser depois.
    h, n, num_incognitas = calcular_h_num(x)
    matriz, vetor = montar_sistema(y, h, num_incognitas)
    matriz, vetor, g = eliminacao_gauss(matriz, vetor, num_incognitas)
    a, b, c, d = calcula_coeficientes(y, h, g, n)

    return a, b, c, d

# ============================================================
# GEOMETRIA DOS PAINÉIS (a partir de uma malha 3D já pronta)
# ============================================================

def formar_paineis(malha):
    num_balizas = len(malha)
    num_linhas_dagua = len(malha[0])
    paineis = []

    for i in range(num_balizas - 1):
        for j in range(num_linhas_dagua - 1):
            p1 = malha[i][j+1]
            p2 = malha[i+1][j+1]
            p3 = malha[i+1][j]
            p4 = malha[i][j]

            if None in (p1, p2, p3, p4):
                continue   # algum canto não existe -> não há painel aqui

            paineis.append([p1, p2, p3, p4])
    return paineis

def vetor_normal_painel(paineis):
    # Calcula o vetor normal de TODOS os painéis de uma vez, via produto
    # vetorial das duas diagonais de cada painel (Eq. 4.3 do relatório).
    paineis_vetor = []
    for i in range(len(paineis)):
        p1 = paineis[i][0]
        p2 = paineis[i][1]
        p3 = paineis[i][2]
        p4 = paineis[i][3]
        diagonal_1 = np.array(p3) - np.array(p1)   # uma diagonal do painel (canto 3 menos canto 1)
        diagonal_2 = np.array(p4) - np.array(p2)   # a outra diagonal (canto 4 menos canto 2)
        normal = np.cross(diagonal_1, diagonal_2)  # produto vetorial das diagonais = vetor normal ao painel
        paineis_vetor.append(normal)
    return paineis_vetor

# ============================================================
# PASSADA VERTICAL: gera pontos dentro de CADA baliza (curva z x y)
# ============================================================

def construir_linha_baliza(baliza, zs_desejados):
    x, z_conhecidos, y_conhecidos = baliza
    pontos_da_baliza = []

    if len(z_conhecidos) == 1:
        z_unico, y_unico = z_conhecidos[0], y_conhecidos[0]
        for z_alvo in zs_desejados:
            if z_alvo == z_unico:
                pontos_da_baliza.append([x, y_unico, z_alvo])
            else:
                pontos_da_baliza.append(None)
        return pontos_da_baliza

    trechos_constantes = encontrar_trechos_constantes(z_conhecidos, y_conhecidos)
    a, b, c, d = ajustar_spline(z_conhecidos, y_conhecidos)

    for z_alvo in zs_desejados:
        if not (z_conhecidos[0] <= z_alvo <= z_conhecidos[-1]):
            pontos_da_baliza.append(None)
            continue
        y_alvo = None
        for z_ini, z_fim, y_const in trechos_constantes:
            if z_ini <= z_alvo <= z_fim:
                y_alvo = y_const
                break
        if y_alvo is None:
            y_alvo = avalia_spline(z_alvo, z_conhecidos, y_conhecidos, a, b, c, d)
        pontos_da_baliza.append([x, y_alvo, z_alvo])

    return pontos_da_baliza

def encontrar_trechos_constantes(x_conhecidos, y_conhecidos, tolerancia=1e-9):
    # Identifica trechos onde 3 ou mais pontos consecutivos compartilham o
    # mesmo valor de y: nesses trechos, o dado real já garante que a forma
    # é reta, então usamos o valor conhecido diretamente, em vez de avaliar
    # a spline (que pode "vazar" curvatura das transições vizinhas devido
    # ao sistema global acoplado).
    trechos = []
    n = len(y_conhecidos)
    i = 0
    while i < n - 1:
        j = i
        while j < n - 1 and abs(y_conhecidos[j+1] - y_conhecidos[i]) < tolerancia:
            j += 1
        if j - i >= 2:   # pelo menos 3 pontos iguais consecutivos
            trechos.append((x_conhecidos[i], x_conhecidos[j], y_conhecidos[i]))
        i = j if j > i else i + 1
    return trechos

# ============================================================
# PASSADA HORIZONTAL: gera balizas novas em CADA altura (curva x x y)
# ============================================================

def construir_linha_dagua(malha_balizas, j, xs_desejados):
    xs_conhecidos, ys_conhecidos = [], []
    z_desta_altura = None
    for pontos_baliza in malha_balizas:
        ponto = pontos_baliza[j]
        if ponto is not None:
            x_i, y_i, z_i = ponto
            xs_conhecidos.append(x_i)
            ys_conhecidos.append(y_i)
            z_desta_altura = z_i
    if len(xs_conhecidos) < 2:
        return [None] * len(xs_desejados)

    trechos_constantes = encontrar_trechos_constantes(xs_conhecidos, ys_conhecidos)
    a, b, c, d = ajustar_spline(xs_conhecidos, ys_conhecidos)

    pontos_da_linha = []
    for x_alvo in xs_desejados:
        if not (xs_conhecidos[0] <= x_alvo <= xs_conhecidos[-1]):
            pontos_da_linha.append(None)
            continue
        y_alvo = None
        for x_ini, x_fim, y_const in trechos_constantes:
            if x_ini <= x_alvo <= x_fim:
                y_alvo = y_const
                break
        if y_alvo is None:
            y_alvo = avalia_spline(x_alvo, xs_conhecidos, ys_conhecidos, a, b, c, d)
        pontos_da_linha.append([x_alvo, y_alvo, z_desta_altura])
    return pontos_da_linha

def gerar_grid_refinado(valores_reais, subdivisoes_por_intervalo=4):
    # Monta uma lista de pontos que passa exatamente pelos valores reais
    # fornecidos e subdivide cada intervalo entre eles, preservando a
    # resolução relativa que os próprios dados já tinham: intervalos reais
    # pequenos (onde a forma muda rápido) continuam gerando intervalos
    # pequenos; intervalos grandes continuam gerando intervalos grandes.
    valores_ordenados = sorted(set(valores_reais))
    grid_refinado = [valores_ordenados[0]]

    for i in range(len(valores_ordenados) - 1):
        inicio = valores_ordenados[i]
        fim = valores_ordenados[i + 1]
        pontos_intervalo = np.linspace(inicio, fim, subdivisoes_por_intervalo + 1)
        grid_refinado.extend(pontos_intervalo[1:].tolist())  # não repete "inicio", já está na lista

    return grid_refinado

def fechar_fundo(malha_final):
    # Gera os painéis que fecham o fundo do casco na quilha (z=0),
    # ligando o ponto de estibordo (y positivo) ao ponto espelhado de
    # bombordo (y negativo) da mesma altura, cobrindo a parte plana da
    # linha d'água zero descrita na seção 3.5 do relatório.
    num_balizas = len(malha_final)
    paineis_fundo = []

    for i in range(num_balizas - 1):
        p_estibordo_i = malha_final[i][0]      # ponto na quilha (z=0), baliza i
        p_estibordo_i1 = malha_final[i + 1][0]  # ponto na quilha, baliza i+1

        if p_estibordo_i is None or p_estibordo_i1 is None:
            continue   # falta dado de quilha numa das balizas -> não há painel aqui

        p_bombordo_i = [p_estibordo_i[0], -p_estibordo_i[1], p_estibordo_i[2]]
        p_bombordo_i1 = [p_estibordo_i1[0], -p_estibordo_i1[1], p_estibordo_i1[2]]

        painel = [p_estibordo_i, p_estibordo_i1, p_bombordo_i1, p_bombordo_i]
        paineis_fundo.append(painel)

    return paineis_fundo

def espelhar_paineis(paineis):
    # Gera o bordo espelhado (bombordo) do casco: para cada painel do
    # bordo já calculado (boreste, y positivo), cria um painel simétrico
    # com y invertido. Troca os pontos 2 e 4 de posição (índices 1 e 3
    # da lista, já que a indexação começa em 0) para manter o vetor
    # normal apontando para fora do casco, conforme a Fig. 3.4.
    paineis_espelhados = []
    for p1, p2, p3, p4 in paineis:
        p1_esp = [p1[0], -p1[1], p1[2]]
        p2_esp = [p2[0], -p2[1], p2[2]]
        p3_esp = [p3[0], -p3[1], p3[2]]
        p4_esp = [p4[0], -p4[1], p4[2]]
        paineis_espelhados.append([p1_esp, p4_esp, p3_esp, p2_esp])  # 2 e 4 trocados
    return paineis_espelhados

def indice_coluna_valida(malha_final, a_partir_do_fim, minimo_pontos=2):
    # Procura, a partir de uma ponta da malha, a primeira coluna com
    # dado suficiente (pelo menos minimo_pontos reais) para formar um
    # fechamento. Isso permite que o código funcione tanto quando a
    # ponta tem dado completo (afinando naturalmente até quase um
    # ponto) quanto quando falta dado (aí ele recua para uma coluna
    # anterior que tenha informação real).
    indices = range(len(malha_final) - 1, -1, -1) if a_partir_do_fim else range(len(malha_final))
    for i in indices:
        pontos_reais = sum(1 for p in malha_final[i] if p is not None)
        if pontos_reais >= minimo_pontos:
            return i
    return None

def fechar_pontas(malha_final, minimo_pontos=2):
    def fechar_coluna(coluna, inverter):
        paineis = []
        for j in range(len(coluna) - 1):
            p_est_j, p_est_j1 = coluna[j], coluna[j + 1]
            if p_est_j is None or p_est_j1 is None:
                continue
            p_bomb_j = [p_est_j[0], -p_est_j[1], p_est_j[2]]
            p_bomb_j1 = [p_est_j1[0], -p_est_j1[1], p_est_j1[2]]
            if not inverter:
                paineis.append([p_est_j, p_bomb_j, p_bomb_j1, p_est_j1])
            else:
                paineis.append([p_est_j1, p_bomb_j1, p_bomb_j, p_est_j])
        return paineis

    idx_proa = indice_coluna_valida(malha_final, a_partir_do_fim=False, minimo_pontos=minimo_pontos)
    idx_popa = indice_coluna_valida(malha_final, a_partir_do_fim=True, minimo_pontos=minimo_pontos)

    paineis_pontas = []
    if idx_proa is not None:
        paineis_pontas += fechar_coluna(malha_final[idx_proa], inverter=False)
    if idx_popa is not None and idx_popa != idx_proa:
        paineis_pontas += fechar_coluna(malha_final[idx_popa], inverter=True)

    return paineis_pontas, idx_proa, idx_popa

# ============================================================
# MONTAGEM FINAL: transpõe a passada horizontal para o formato malha[i][j]
# ============================================================

def montar_malha_final(malha_linhas_dagua, num_balizas_novas):
    malha_final = []
    for i in range(num_balizas_novas):
        coluna = [linha[i] for linha in malha_linhas_dagua]
        malha_final.append(coluna)
    return malha_final

'''
# ============================================================
# SCRIPT PRINCIPAL: monta a malha de exemplo de ponta a ponta e testa
# ============================================================

pontos_desejados = np.linspace(0, 2, 9)   # sobra de um teste anterior (exercício do professor)
malha_fina = pontos_desejados.tolist()     # não é mais usada no fluxo atual - pode ser removida depois

# --- descobre a faixa global de alturas (para a passada vertical) ---
alturas_maximas = []
for baliza in malha_teste:
    alturas_maximas.append(baliza[1][-1])   # baliza[1] é a lista z_conhecidos; [-1] é o maior valor dela

z_max_global = max(alturas_maximas)          # maior altura entre todas as balizas (aqui: a meia-nau, z=4)
z_min_global = 0.0                           # geralmente a quilha, comum a todas as balizas

zs_desejados = np.linspace(z_min_global, z_max_global, 9).tolist()   # 9 alturas igualmente espaçadas, 0 a 4

# --- PASSADA VERTICAL: roda a spline dentro de cada baliza ---
malha_balizas = []
for baliza in malha_teste:
    pontos = construir_linha_baliza(baliza, zs_desejados)
    malha_balizas.append(pontos)

# --- descobre a faixa global de posições x (para a passada horizontal) ---
xs_balizas = []
for baliza in malha_teste:
    xs_balizas.append(baliza[0])   # baliza[0] é a posição x dela (popa=0, meia-nau=5, proa=10)

x_min_global = min(xs_balizas)
x_max_global = max(xs_balizas)

xs_desejados = np.linspace(x_min_global, x_max_global, 9).tolist()   # 9 posições de baliza, 0 a 10

# --- PASSADA HORIZONTAL: roda a spline entre as balizas, altura por altura ---
malha_linhas_dagua = []
for j in range(len(zs_desejados)):                              # uma chamada para cada altura desejada
    linha = construir_linha_dagua(malha_balizas, j, xs_desejados)
    malha_linhas_dagua.append(linha)

# --- junta as duas passadas no formato que formar_paineis espera ---
malha_final = montar_malha_final(malha_linhas_dagua, len(xs_desejados))

# --- gera os painéis e confere a orientação deles ---
paineis = formar_paineis(malha_final)
paineis_vetor = vetor_normal_painel(paineis)

print([len(coluna) for coluna in malha_final])          # tamanho de cada baliza da malha final (esperado: 9x "5")
print([v.tolist()[1] for v in paineis_vetor])            # só o componente y de cada vetor normal (esperado: só +)
'''