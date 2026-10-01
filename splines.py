import matplotlib.pyplot as plt   # biblioteca de gráficos (importada, ainda não usada no código)
import numpy as np                # usado para np.linspace, np.array e np.cross

TOL_Y = 1e-3

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
    # Forma um painel (quadrilátero) para cada "célula" da malha.
    # Como a malha agora segue o contorno real do casco (não há mais pontos
    # "de enchimento" fora dele), basta pular células com canto None e as
    # que ficariam inteiras em cima da linha de centro (área nula).
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

            if all(abs(p[1]) < TOL_Y for p in (p1, p2, p3, p4)):
                continue   # os 4 cantos na linha de centro -> painel sem área

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

    # --- IDENTIFICAÇÃO DO LIMITE INFERIOR REAL ---
    idx_inicio_real = 0
    for i, y_val in enumerate(y_conhecidos):
        if abs(y_val) > 1e-4:
            idx_inicio_real = max(0, i - 1)
            break

    z_min_real = z_conhecidos[idx_inicio_real]
    z_max_real = z_conhecidos[-1]

    trechos_constantes = encontrar_trechos_constantes(z_conhecidos, y_conhecidos)
    a, b, c, d = ajustar_spline(z_conhecidos, y_conhecidos)

    for z_alvo in zs_desejados:
        # Força o retorno de None para as regiões vazias fora do aço real
        if not (z_min_real <= z_alvo <= z_max_real):
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
        # Se a altura desejada estiver fora do alcance daquela baliza,
        # ela zera (encosta na linha de centro do navio) em vez de dar None.
        # Isso faz o casco afinar nas pontas sem abrir buracos na malha!
        if not (z_conhecidos[0] <= z_alvo <= z_conhecidos[-1]):
            pontos_da_baliza.append([x, 0.0, z_alvo])
            continue
            
        y_alvo = None
        for z_ini, z_fim, y_const in trechos_constantes:
            if z_ini <= z_alvo <= z_fim:
                y_alvo = y_const
                break
        if y_alvo is None:
            y_alvo = avalia_spline(z_alvo, z_conhecidos, y_conhecidos, a, b, c, d)
            
        # Garante que larguras calculadas como negativas fiquem na linha de centro
        if y_alvo < 0:
            y_alvo = 0.0
            
        pontos_da_baliza.append([x, y_alvo, z_alvo])

    return pontos_da_baliza


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
        return [[x, 0.0, 0.0] for x in xs_desejados]

    trechos_constantes = encontrar_trechos_constantes(xs_conhecidos, ys_conhecidos)
    a, b, c, d = ajustar_spline(xs_conhecidos, ys_conhecidos)

    pontos_da_linha = []
    for x_alvo in xs_desejados:
        if not (xs_conhecidos[0] <= x_alvo <= xs_conhecidos[-1]):
            pontos_da_linha.append([x_alvo, 0.0, z_desta_altura])
            continue
        
        y_alvo = None
        for x_ini, x_fim, y_const in trechos_constantes:
            if x_ini <= x_alvo <= x_fim:
                y_alvo = y_const
                break
        if y_alvo is None:
            y_alvo = avalia_spline(x_alvo, xs_conhecidos, ys_conhecidos, a, b, c, d)
            
        if y_alvo < 0:
            y_alvo = 0.0
            
        pontos_da_linha.append([x_alvo, y_alvo, z_desta_altura])
            
    return pontos_da_linha

def estimar_ponta(x_ext, y_ext, x_int, y_int, x_zero, sentido):
    # Estima onde a linha d'água chega a y = 0, a partir do último ponto com
    # y > 0 (x_ext, y_ext) e do vizinho interno (x_int, y_int): prolonga a
    # inclinação local até y = 0, sem passar da estação vizinha onde a tabela
    # já mostra y = 0 (x_zero). sentido = +1 rumo ao maior x, -1 rumo ao menor x.
    limite = abs(x_zero - x_ext)
    if x_int is None:
        return x_ext + sentido * limite
    inclinacao = (y_ext - y_int) / ((x_ext - x_int) * sentido)   # dy/dx indo para fora
    if inclinacao >= 0:
        return x_ext + sentido * limite
    distancia = min(y_ext / (-inclinacao), limite)
    return x_ext + sentido * distancia

def posicao_ponta_por_z(lista_ext, z):
    # lista_ext = [(z, x_ponta), ...] em ordem crescente de z: posição (x) onde a
    # linha d'água da altura z fecha na linha de centro (contorno de proa ou popa,
    # lido do plano de linhas). Interpola linearmente; fora da faixa devolve None.
    if not lista_ext or z < lista_ext[0][0] - 1e-9 or z > lista_ext[-1][0] + 1e-9:
        return None
    for (z0, x0), (z1, x1) in zip(lista_ext, lista_ext[1:]):
        if z0 - 1e-9 <= z <= z1 + 1e-9:
            return x0 + (x1 - x0) * (z - z0) / (z1 - z0) if z1 > z0 else x0
    return lista_ext[-1][1]

def construir_linha_dagua_contorno(malha_balizas, j, n_colunas, concentracao_pontas=0.0,
                                   ext_ini=None, ext_fim=None):
    # Gera os n_colunas pontos da linha d'água j, entre as suas duas pontas.
    # concentracao_pontas (0 a 1): 0 = pontos uniformes; valores maiores
    # adensam os pontos perto das extremidades, onde a curvatura é maior.
    # ext_ini / ext_fim (opcionais): contorno além da primeira / última estação da
    # tabela, como lista [(z, x_ponta), ...]. Só é usado nas linhas d'água cujo
    # casco continua além dessa estação (y > 0 nela), como o "overhang" da proa
    # acima do calado e da popa de cruzador.
    xs, ys, z_altura = [], [], None
    for pontos_baliza in malha_balizas:
        ponto = pontos_baliza[j]
        if ponto is not None:
            xs.append(ponto[0])
            ys.append(ponto[1])
            z_altura = ponto[2]

    positivos = [i for i, y in enumerate(ys) if y > TOL_Y]
    if len(xs) < 2 or not positivos:
        return None   # esta altura não corta o casco

    i_a, i_b = positivos[0], positivos[-1]     # primeira e última estação com casco
    nos_x, nos_y = xs[i_a:i_b + 1], ys[i_a:i_b + 1]

    # pontas: se existe uma estação vizinha com y = 0, a linha d'água fecha ali;
    # senão o casco continua além da tabela: usa o contorno informado (se houver)
    # ou deixa a ponta "aberta" (parede plana na última estação).
    ponta_ini = ponta_fim = None
    if i_a > 0:
        x_int = nos_x[1] if len(nos_x) > 1 else None
        y_int = nos_y[1] if len(nos_x) > 1 else None
        ponta_ini = estimar_ponta(nos_x[0], nos_y[0], x_int, y_int, xs[i_a - 1], -1)
    else:
        ponta_ini = posicao_ponta_por_z(ext_ini, z_altura)
        if ponta_ini is not None and ponta_ini >= nos_x[0] - 1e-6:
            ponta_ini = None
    if i_b < len(xs) - 1:
        x_int = nos_x[-2] if len(nos_x) > 1 else None
        y_int = nos_y[-2] if len(nos_x) > 1 else None
        ponta_fim = estimar_ponta(nos_x[-1], nos_y[-1], x_int, y_int, xs[i_b + 1], +1)
    else:
        ponta_fim = posicao_ponta_por_z(ext_fim, z_altura)
        if ponta_fim is not None and ponta_fim <= nos_x[-1] + 1e-6:
            ponta_fim = None

    fecha_ini = ponta_ini is not None and abs(ponta_ini - nos_x[0]) > 1e-6
    fecha_fim = ponta_fim is not None and abs(ponta_fim - nos_x[-1]) > 1e-6
    if fecha_ini:
        nos_x.insert(0, ponta_ini); nos_y.insert(0, 0.0)
    if fecha_fim:
        nos_x.append(ponta_fim); nos_y.append(0.0)

    a, b, c, d = ajustar_spline(nos_x, nos_y)
    trechos_constantes = encontrar_trechos_constantes(nos_x, nos_y)
    x_ini, x_fim = nos_x[0], nos_x[-1]

    pontos_da_linha = []
    for k in range(n_colunas):
        u = k / (n_colunas - 1)
        w = (1 - concentracao_pontas) * u + concentracao_pontas * (0.5 - 0.5 * np.cos(np.pi * u))
        x_alvo = x_ini + (x_fim - x_ini) * w
        if k == 0:
            x_alvo = x_ini                 # evita erro de arredondamento nas pontas
        if k == n_colunas - 1:
            x_alvo = x_fim

        y_alvo = None
        for x_i, x_f, y_const in trechos_constantes:
            if x_i <= x_alvo <= x_f:
                y_alvo = y_const
                break
        if y_alvo is None:
            y_alvo = avalia_spline(x_alvo, nos_x, nos_y, a, b, c, d)
        if y_alvo < 0:
            y_alvo = 0.0
        if (k == 0 and fecha_ini) or (k == n_colunas - 1 and fecha_fim):
            y_alvo = 0.0                   # a ponta fecha exatamente na linha de centro
        pontos_da_linha.append([x_alvo, y_alvo, z_altura])
    return pontos_da_linha

def construir_malha_por_contorno(malha_balizas, n_colunas, concentracao_pontas=0.0,
                                 ext_ini=None, ext_fim=None):
    # Monta a malha final malha[i][j] (i = coluna ao longo do comprimento,
    # j = linha d'água) já no contorno real. ext_ini / ext_fim: ver
    # construir_linha_dagua_contorno. Substitui o trio
    # construir_linha_dagua + montar_malha_final + aplicar_mascara_serie60.
    linhas = []
    for j in range(len(malha_balizas[0])):
        linha = construir_linha_dagua_contorno(malha_balizas, j, n_colunas, concentracao_pontas,
                                               ext_ini, ext_fim)
        if linha is None:
            print(f"Aviso: a linha d'água de índice {j} não corta o casco e foi ignorada.")
            continue
        linhas.append(linha)
    return montar_malha_final(linhas, n_colunas)

def ponto_quilha(p):
    # Ponto da linha de centro que fecha o fundo abaixo do ponto p (a primeira
    # linha d'água de uma coluna). Se p já está na linha de centro (ponta do
    # casco), não há o que fechar e devolve o próprio ponto.
    if abs(p[1]) < TOL_Y:
        return [p[0], 0.0, p[2]]
    return [p[0], 0.0, 0.0]

def contorno_perfil(malha_final):
    # Silhueta do casco visto de lado (plano de simetria, y = 0): quilha,
    # contorno da extremidade final, convés e contorno da extremidade inicial.
    quilha = [ponto_quilha(col[0]) for col in malha_final]
    fim = [[p[0], 0.0, p[2]] for p in malha_final[-1]]
    conves = [[col[-1][0], 0.0, col[-1][2]] for col in malha_final]
    ini = [[p[0], 0.0, p[2]] for p in malha_final[0]]
    contorno = quilha + fim + conves[::-1] + ini[::-1] + [quilha[0]]
    return contorno

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
    # Fecha o fundo: liga a primeira linha d'água de cada coluna (bordo de
    # boreste e de bombordo) a uma quilha central na linha de centro.
    # Nas pontas (onde o ponto já está na linha de centro) o painel vira um
    # triângulo ou some, sem criar "aletas" artificiais.
    num_balizas = len(malha_final)
    paineis_fundo = []

    for i in range(num_balizas - 1):
        p_estibordo_i = malha_final[i][0]
        p_estibordo_i1 = malha_final[i + 1][0]

        if p_estibordo_i is None or p_estibordo_i1 is None:
            continue
        if abs(p_estibordo_i[1]) < TOL_Y and abs(p_estibordo_i1[1]) < TOL_Y:
            continue   # fundo sem largura entre as duas colunas

        p_quilha_i = ponto_quilha(p_estibordo_i)
        p_quilha_i1 = ponto_quilha(p_estibordo_i1)

        paineis_fundo.append([p_estibordo_i, p_estibordo_i1, p_quilha_i1, p_quilha_i])

        p_bombordo_i = [p_estibordo_i[0], -p_estibordo_i[1], p_estibordo_i[2]]
        p_bombordo_i1 = [p_estibordo_i1[0], -p_estibordo_i1[1], p_estibordo_i1[2]]
        paineis_fundo.append([p_quilha_i, p_quilha_i1, p_bombordo_i1, p_bombordo_i])

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
    def fechar_coluna_lamina(coluna, inverter):
        paineis = []
        for j in range(len(coluna) - 1):
            p_est_j, p_est_j1 = coluna[j], coluna[j + 1]
            if p_est_j is None or p_est_j1 is None:
                continue
            
            x_j, y_j, z_j = p_est_j
            x_j1, y_j1, z_j1 = p_est_j1

            p_centro_j = [x_j, 0.0, z_j]
            p_centro_j1 = [x_j1, 0.0, z_j1]
            
            # Se a própria lateral do navio já convergiu para zero na tabela,
            # significa que a rampa naval fechou sozinha. Não cria painel extra!
            if abs(y_j) < 1e-3 and abs(y_j1) < 1e-3:
                continue

            p_bomb_j = [x_j, -y_j, z_j]
            p_bomb_j1 = [x_j1, -y_j1, z_j1]

            if not inverter:
                paineis.append([p_est_j, p_centro_j, p_centro_j1, p_est_j1])
                paineis.append([p_centro_j, p_bomb_j, p_bomb_j1, p_centro_j1])
            else:
                paineis.append([p_est_j1, p_centro_j1, p_centro_j, p_est_j])
                paineis.append([p_bomb_j1, p_centro_j1, p_centro_j, p_bomb_j])
        return paineis

    # Forçamos o fechamento sempre nos limites físicos da malha gerada
    paineis_pontas = []
    paineis_pontas += fechar_coluna_lamina(malha_final[0], inverter=False)
    paineis_pontas += fechar_coluna_lamina(malha_final[-1], inverter=True)

    return paineis_pontas, 0, len(malha_final) - 1

def montar_malha_final(malha_linhas_dagua, num_balizas_novas):
    malha_final = []
    for i in range(num_balizas_novas):
        coluna = [linha[i] for linha in malha_linhas_dagua]
        malha_final.append(coluna)
    return malha_final
