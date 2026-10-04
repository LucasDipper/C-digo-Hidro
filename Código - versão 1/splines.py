import numpy as np                # Biblioteca numérica: usada aqui em np.array, np.cross e np.linspace

# ============================================================
# VISÃO GERAL E CONVENÇÕES DO ARQUIVO
# Este arquivo transforma as balizas medidas (tabela de cotas) em uma malha 3D suave
# e depois em painéis quadriláteros que formam a superfície do casco.
# Fluxo geral:
#   1) passada vertical   -> spline y(z) dentro de cada baliza
#   2) passada horizontal -> spline y(x) em cada altura z
#   3) montagem da malha  -> malha_final[i][j]
#   4) painéis            -> boreste, bombordo (espelhado), fundo e pontas
# Convenções usadas em todo o arquivo:
#   - um ponto é uma lista [x, y, z]; "None" significa "sem dado naquela posição"
#   - malha[i][j]: i = baliza (direção x; 0 = proa, último = popa), j = linha d'água (direção z)
#   - y > 0 é o lado de boreste; bombordo é obtido espelhando (y -> -y)
#   - os vetores x e z conhecidos são assumidos em ordem crescente
# ============================================================

# ============================================================
# MOTOR DE SPLINE CÚBICA NATURAL (curva única: um x, um y)
# Resolve a spline de uma curva por vez. É chamado por ajustar_spline(), que roda o
# processo completo: espaçamentos -> sistema linear -> Gauss -> coeficientes.
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
    # Cada linha do sistema é a condição de continuidade da primeira derivada no nó k.
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
        # (diferença entre as inclinações dos dois trechos vizinhos ao nó k)
        vetor[linha] = 6 * ((y[k + 1] - y[k]) / h[k] - (y[k] - y[k - 1]) / h[k - 1])

    return matriz, vetor   # sistema pronto para ser resolvido pela eliminação de Gauss

def eliminacao_gauss(matriz, vetor, num_incognitas):
    # Resolve o sistema linear (matriz, vetor) por eliminação de Gauss:
    # primeiro escalona (zera abaixo da diagonal), depois resolve de baixo para cima.
    # Não usa pivotamento: não é necessário aqui, pois a matriz é tridiagonal e
    # diagonalmente dominante (o pivô nunca é zero).

    # --- ESCALONAMENTO ---
    for coluna in range(num_incognitas - 1):                # percorre cada coluna a ser zerada
        for l_alvo in range(coluna + 1, num_incognitas):     # percorre as linhas abaixo da diagonal nessa coluna
            fator = matriz[l_alvo][coluna] / matriz[coluna][coluna]  # fator que zera o elemento
            vetor[l_alvo] -= fator * vetor[coluna]            # aplica a mesma operação no vetor

            for i in range(coluna, num_incognitas):           # aplica a subtração na linha da matriz, da coluna atual até o fim
                matriz[l_alvo][i] -= fator * matriz[coluna][i]

    # --- SUBSTITUIÇÃO REVERSA ---
    g = [None] * num_incognitas                                       # vai guardar a solução (os g_k)
    for i in range(num_incognitas - 1, -1, -1):  # percorre as linhas de baixo para cima: do último índice até 0
        soma = 0                                  # acumula a contribuição das incógnitas já resolvidas
        for j in range(num_incognitas - 1, -1, -1):
            if j > i:                              # só entram as incógnitas à direita de i (já resolvidas antes)
                soma += g[j]*matriz[i][j]
        g[i] = (vetor[i] - soma)/matriz[i][i]       # isola g[i] na equação

    return matriz, vetor, g   # matriz/vetor já escalonados (não usados depois) + a solução g

def calcula_coeficientes(y, h, g, n):
    # Calcula os coeficientes a, b, c, d de cada trecho da spline, a partir dos g_k já resolvidos.
    # Em cada trecho k (entre x[k-1] e x[k]) o polinômio é escrito em potências de (x - x[k]):
    # S_k(x) = a*(x-x_k)^3 + b*(x-x_k)^2 + c*(x-x_k) + d
    g_completo = [0] + g + [0]   # reinsere g0=0 e g(n-1)=0, removidos do sistema pela condição natural
    a = [None] * n               # coeficiente do termo cúbico de cada trecho
    b = [None] * n               # coeficiente do termo quadrático
    c = [None] * n               # coeficiente do termo linear
    d = [None] * n               # termo constante (valor no nó da direita do trecho)
    # as listas têm tamanho n, mas a posição 0 fica vazia (None): o trecho k vai de x[k-1] a x[k], com k >= 1

    for k in range(1, n):        # um conjunto de coeficientes para cada trecho (k = 1 até n-1)
        a[k] = (g_completo[k] - g_completo[k-1])/(6*h[k-1])   # variação da 2ª derivada ao longo do trecho
        b[k] = g_completo[k]/2                                  # metade da 2ª derivada no nó da direita
        c[k] = ((y[k] - y[k-1])/h[k-1]) + ((2*h[k-1]*g_completo[k] + g_completo[k-1]*h[k-1])/6)  # inclinação no nó da direita
        d[k] = y[k]               # o termo constante é sempre o valor de y no nó da direita

    return a, b, c, d

def avalia_spline(x_alvo, x, y, a, b, c, d):
    # Avalia a spline já ajustada num ponto x_alvo qualquer: primeiro acha
    # em qual trecho [x[k-1], x[k]] ele cai, depois calcula o polinômio ali.
    # Obs.: se x_alvo estiver fora do intervalo dos nós, o laço termina em k = n-1 e a função
    # extrapola com o último trecho; por isso quem chama confere o intervalo antes.
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
    h, n, num_incognitas = calcular_h_num(x)                          # espaçamentos e tamanhos
    matriz, vetor = montar_sistema(y, h, num_incognitas)              # sistema A*g = vetor
    matriz, vetor, g = eliminacao_gauss(matriz, vetor, num_incognitas) # resolve o sistema -> g_k
    a, b, c, d = calcula_coeficientes(y, h, g, n)                     # coeficientes de cada trecho

    return a, b, c, d

# ============================================================
# GEOMETRIA DOS PAINÉIS (a partir de uma malha 3D já pronta)
# Aqui nascem os painéis do lado de boreste e os vetores normais de cada painel.
# O espelhamento (bombordo) e os fechamentos (fundo e pontas) estão mais abaixo no arquivo.
# ============================================================

def formar_paineis(malha):
    # Forma um painel quadrilátero para cada "quadradinho" da malha, entre duas balizas vizinhas
    # (i e i+1) e duas linhas d'água vizinhas (j e j+1). Gera só o lado de boreste (y > 0);
    # o lado de bombordo vem de espelhar_paineis().
    num_balizas = len(malha)            # nº de balizas (direção x)
    num_linhas_dagua = len(malha[0])    # nº de linhas d'água (direção z), usando a 1ª baliza como referência
    paineis = []                        # lista final; cada painel = [p1, p2, p3, p4]

    for i in range(num_balizas - 1):            # n balizas -> n-1 faixas entre elas
        for j in range(num_linhas_dagua - 1):   # n linhas d'água -> n-1 faixas entre elas
            # quatro cantos em sequência ao redor do painel (a ordem define o sentido da normal):
            p1 = malha[i][j+1]      # baliza i,   linha d'água de cima
            p2 = malha[i+1][j+1]    # baliza i+1, linha d'água de cima
            p3 = malha[i+1][j]      # baliza i+1, linha d'água de baixo
            p4 = malha[i][j]        # baliza i,   linha d'água de baixo

            if None in (p1, p2, p3, p4):
                continue   # algum canto não existe -> não há painel aqui

            paineis.append([p1, p2, p3, p4])   # guarda o painel completo
    return paineis

def vetor_normal_painel(paineis):
    # Calcula o vetor normal de TODOS os painéis de uma vez, via produto
    # vetorial das duas diagonais de cada painel.
    # O vetor NÃO é normalizado: seu módulo é igual a 2x a área do painel.
    paineis_vetor = []
    for i in range(len(paineis)):
        p1 = paineis[i][0]   # os quatro cantos do painel i
        p2 = paineis[i][1]
        p3 = paineis[i][2]
        p4 = paineis[i][3]
        diagonal_1 = np.array(p3) - np.array(p1)   # uma diagonal do painel (canto 3 menos canto 1)
        diagonal_2 = np.array(p4) - np.array(p2)   # a outra diagonal (canto 4 menos canto 2)
        normal = np.cross(diagonal_1, diagonal_2)  # produto vetorial das diagonais = vetor normal ao painel
        paineis_vetor.append(normal)               # guarda na mesma ordem dos painéis
    return paineis_vetor

# ============================================================
# PASSADA VERTICAL: gera pontos dentro de CADA baliza (curva y em função de z)
# Para cada baliza original, ajusta uma spline y(z) e avalia nas alturas do grid refinado.
# ============================================================

def construir_linha_baliza(baliza, zs_desejados):
    # Recebe UMA baliza no formato (x, lista de z, lista de y) e devolve a lista de pontos
    # [x, y, z] dessa baliza em cada altura de zs_desejados (None onde não há dado).
    x, z_conhecidos, y_conhecidos = baliza   # desempacota: posição x, alturas medidas, cotas (y) medidas
    pontos_da_baliza = []                    # resultado: um ponto (ou None) por altura desejada

    # Caso especial: baliza com um único ponto medido -> não dá para ajustar spline.
    # Só existe ponto na própria altura medida; nas demais alturas fica None.
    # (a comparação z_alvo == z_unico funciona porque o grid refinado contém os valores reais exatos)
    if len(z_conhecidos) == 1:
        z_unico, y_unico = z_conhecidos[0], y_conhecidos[0]
        for z_alvo in zs_desejados:
            if z_alvo == z_unico:
                pontos_da_baliza.append([x, y_unico, z_alvo])
            else:
                pontos_da_baliza.append(None)
        return pontos_da_baliza

    # Caso geral: baliza com 2 ou mais pontos medidos
    trechos_constantes = encontrar_trechos_constantes(z_conhecidos, y_conhecidos)  # trechos retos já conhecidos
    a, b, c, d = ajustar_spline(z_conhecidos, y_conhecidos)                         # spline y(z) desta baliza

    for z_alvo in zs_desejados:
        if not (z_conhecidos[0] <= z_alvo <= z_conhecidos[-1]):
            pontos_da_baliza.append(None)   # altura fora do intervalo medido: não extrapola
            continue
        y_alvo = None
        for z_ini, z_fim, y_const in trechos_constantes:
            if z_ini <= z_alvo <= z_fim:
                y_alvo = y_const            # dentro de um trecho constante: usa o valor real direto
                break
        if y_alvo is None:
            # fora dos trechos constantes: avalia a spline
            y_alvo = avalia_spline(z_alvo, z_conhecidos, y_conhecidos, a, b, c, d)
        pontos_da_baliza.append([x, y_alvo, z_alvo])

    return pontos_da_baliza

# ============================================================
# AUXILIAR (usada pelas duas passadas): detecção de trechos constantes
# ============================================================

def encontrar_trechos_constantes(x_conhecidos, y_conhecidos, tolerancia=1e-9):
    # Identifica trechos onde 3 ou mais pontos consecutivos compartilham o
    # mesmo valor de y: nesses trechos, o dado real já garante que a forma
    # é reta, então usamos o valor conhecido diretamente, em vez de avaliar
    # a spline (que pode "vazar" curvatura das transições vizinhas devido
    # ao sistema global acoplado).
    # Serve tanto para a passada vertical (x_conhecidos = alturas z) quanto para a horizontal
    # (x_conhecidos = posições x). Devolve uma lista de tuplas (inicio, fim, y_constante).
    trechos = []
    n = len(y_conhecidos)
    i = 0                                    # i = início do trecho em análise
    while i < n - 1:
        j = i                                # j = último ponto do trecho que começa em i
        # avança j enquanto o próximo y for igual ao y do início do trecho (dentro da tolerância)
        while j < n - 1 and abs(y_conhecidos[j+1] - y_conhecidos[i]) < tolerancia:
            j += 1
        if j - i >= 2:   # pelo menos 3 pontos iguais consecutivos
            trechos.append((x_conhecidos[i], x_conhecidos[j], y_conhecidos[i]))
        i = j if j > i else i + 1   # continua a busca a partir do fim do trecho (ou do próximo ponto, se não houve igualdade)
    return trechos

# ============================================================
# PASSADA HORIZONTAL: gera balizas novas em CADA altura (curva y em função de x)
# Para cada altura z do grid refinado, ajusta uma spline y(x) através das balizas
# originais e avalia nas posições x desejadas.
# ============================================================

def construir_linha_dagua(malha_balizas, j, xs_desejados):
    # Constrói a linha d'água de índice j: lista de pontos [x, y, z] nas posições xs_desejados
    # (None onde não há dado).
    xs_conhecidos, ys_conhecidos = [], []   # pontos reais (x, y) disponíveis nesta altura
    z_desta_altura = None
    for pontos_baliza in malha_balizas:      # percorre todas as balizas
        ponto = pontos_baliza[j]             # ponto da baliza na altura j
        if ponto is not None:                # balizas sem dado nessa altura são ignoradas
            x_i, y_i, z_i = ponto
            xs_conhecidos.append(x_i)
            ys_conhecidos.append(y_i)
            z_desta_altura = z_i             # z é o mesmo em todos os pontos desta altura
    if len(xs_conhecidos) < 2:
        return [None] * len(xs_desejados)    # menos de 2 pontos: não dá para ajustar curva

    trechos_constantes = encontrar_trechos_constantes(xs_conhecidos, ys_conhecidos)  # trechos retos ao longo de x
    a, b, c, d = ajustar_spline(xs_conhecidos, ys_conhecidos)                        # spline y(x) desta altura

    pontos_da_linha = []
    for x_alvo in xs_desejados:
        if not (xs_conhecidos[0] <= x_alvo <= xs_conhecidos[-1]):
            # fora do intervalo com dado: não extrapola (é aqui que surgem os avisos de proa/popa no main)
            pontos_da_linha.append(None)
            continue
        y_alvo = None
        for x_ini, x_fim, y_const in trechos_constantes:
            if x_ini <= x_alvo <= x_fim:
                y_alvo = y_const             # dentro de um trecho constante: usa o valor real direto
                break
        if y_alvo is None:
            # fora dos trechos constantes: avalia a spline
            y_alvo = avalia_spline(x_alvo, xs_conhecidos, ys_conhecidos, a, b, c, d)
        pontos_da_linha.append([x_alvo, y_alvo, z_desta_altura])
    return pontos_da_linha

# ============================================================
# GRID REFINADO: cria as alturas (ou posições) usadas para avaliar as splines
# ============================================================

def gerar_grid_refinado(valores_reais, subdivisoes_por_intervalo=4):
    # Monta uma lista de pontos que passa exatamente pelos valores reais
    # fornecidos e subdivide cada intervalo entre eles, preservando a
    # resolução relativa que os próprios dados já tinham: intervalos reais
    # pequenos (onde a forma muda rápido) continuam gerando intervalos
    # pequenos; intervalos grandes continuam gerando intervalos grandes.
    valores_ordenados = sorted(set(valores_reais))   # remove repetidos e ordena do menor para o maior
    grid_refinado = [valores_ordenados[0]]           # começa pelo primeiro valor real

    for i in range(len(valores_ordenados) - 1):      # percorre cada intervalo entre valores reais vizinhos
        inicio = valores_ordenados[i]
        fim = valores_ordenados[i + 1]
        # divide o intervalo em "subdivisoes_por_intervalo" partes iguais (inclui início e fim)
        pontos_intervalo = np.linspace(inicio, fim, subdivisoes_por_intervalo + 1)
        grid_refinado.extend(pontos_intervalo[1:].tolist())  # não repete "inicio", já está na lista

    return grid_refinado

# ============================================================
# FECHAMENTO E ESPELHAMENTO DOS PAINÉIS
# Complementam os painéis de boreste: fundo (quilha), lado de bombordo e pontas (proa e popa).
# ============================================================

def fechar_fundo(malha_final):
    # Gera os painéis que fecham o fundo do casco na quilha (z=0),
    # ligando o ponto de boreste (y positivo) ao ponto espelhado de
    # bombordo (y negativo) da mesma altura, cobrindo a parte plana da
    # linha d'água zero descrita na seção 3.5 do relatório.
    # Usa a linha d'água de índice 0 (a mais baixa do grid, a quilha). Se y = 0 nessa linha,
    # o painel colapsa em uma linha (área nula).
    num_balizas = len(malha_final)
    paineis_fundo = []

    for i in range(num_balizas - 1):             # cada faixa entre duas balizas vizinhas
        p_estibordo_i = malha_final[i][0]        # ponto na quilha (z=0), baliza i (lado de boreste)
        p_estibordo_i1 = malha_final[i + 1][0]   # ponto na quilha, baliza i+1

        if p_estibordo_i is None or p_estibordo_i1 is None:
            continue   # falta dado de quilha numa das balizas -> não há painel aqui

        # pontos espelhados para bombordo (y com sinal trocado)
        p_bombordo_i = [p_estibordo_i[0], -p_estibordo_i[1], p_estibordo_i[2]]
        p_bombordo_i1 = [p_estibordo_i1[0], -p_estibordo_i1[1], p_estibordo_i1[2]]

        # painel atravessando o casco de boreste a bombordo, nas duas balizas
        painel = [p_estibordo_i, p_estibordo_i1, p_bombordo_i1, p_bombordo_i]
        paineis_fundo.append(painel)

    return paineis_fundo

def espelhar_paineis(paineis):
    # Gera o lado oposto (bombordo) do casco: para cada painel de boreste
    # já calculado (y positivo), cria um painel simétrico com y invertido.
    # Espelhar inverte o sentido de percurso dos cantos; por isso troca os pontos 2 e 4
    # de posição (índices 1 e 3 da lista, já que a indexação começa em 0) para manter
    # o vetor normal apontando para fora do casco, conforme a Fig. 3.4.
    paineis_espelhados = []
    for p1, p2, p3, p4 in paineis:            # desempacota os quatro cantos de cada painel
        p1_esp = [p1[0], -p1[1], p1[2]]       # mesmo x e z, y com sinal trocado
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
    # a_partir_do_fim=True procura da popa (última baliza) para a proa; False, da proa para a popa.
    # Devolve o índice da baliza encontrada, ou None se nenhuma tiver dado suficiente.
    indices = range(len(malha_final) - 1, -1, -1) if a_partir_do_fim else range(len(malha_final))
    for i in indices:
        pontos_reais = sum(1 for p in malha_final[i] if p is not None)   # conta os pontos que existem na coluna i
        if pontos_reais >= minimo_pontos:
            return i
    return None

def fechar_pontas(malha_final, minimo_pontos=2):
    # Fecha a proa e a popa com painéis que atravessam o casco de boreste a bombordo.
    # Devolve os painéis das pontas e os índices das balizas onde o fechamento foi feito
    # (o main usa esses índices para avisar quando faltou dado nas pontas).

    def fechar_coluna(coluna, inverter):
        # Gera os painéis de UMA baliza (coluna de pontos de baixo para cima), ligando cada par de
        # linhas d'água vizinhas ao seu espelho de bombordo.
        # "inverter" troca a ordem dos cantos: proa e popa olham para lados opostos, então a
        # ordem precisa ser oposta para que a normal aponte para fora do casco nas duas pontas.
        paineis = []
        for j in range(len(coluna) - 1):
            p_est_j, p_est_j1 = coluna[j], coluna[j + 1]   # dois pontos de boreste vizinhos (alturas j e j+1)
            if p_est_j is None or p_est_j1 is None:
                continue                                    # falta dado -> sem painel nesse trecho
            p_bomb_j = [p_est_j[0], -p_est_j[1], p_est_j[2]]      # correspondentes espelhados em bombordo
            p_bomb_j1 = [p_est_j1[0], -p_est_j1[1], p_est_j1[2]]
            if not inverter:
                paineis.append([p_est_j, p_bomb_j, p_bomb_j1, p_est_j1])
            else:
                paineis.append([p_est_j1, p_bomb_j1, p_bomb_j, p_est_j])
        return paineis

    # índices das balizas válidas mais próximas da proa (início) e da popa (fim)
    idx_proa = indice_coluna_valida(malha_final, a_partir_do_fim=False, minimo_pontos=minimo_pontos)
    idx_popa = indice_coluna_valida(malha_final, a_partir_do_fim=True, minimo_pontos=minimo_pontos)

    paineis_pontas = []
    if idx_proa is not None:
        paineis_pontas += fechar_coluna(malha_final[idx_proa], inverter=False)   # fechamento da proa
    if idx_popa is not None and idx_popa != idx_proa:                             # evita fechar duas vezes a mesma baliza
        paineis_pontas += fechar_coluna(malha_final[idx_popa], inverter=True)    # fechamento da popa

    return paineis_pontas, idx_proa, idx_popa

# ============================================================
# MONTAGEM FINAL: transpõe a passada horizontal para o formato malha[i][j]
# A passada horizontal entrega uma lista de linhas d'água (índice = altura j);
# aqui ela é reorganizada em uma lista de balizas (índice = baliza i).
# ============================================================

def montar_malha_final(malha_linhas_dagua, num_balizas_novas):
    malha_final = []
    for i in range(num_balizas_novas):                          # uma coluna para cada baliza nova
        coluna = [linha[i] for linha in malha_linhas_dagua]     # i-ésimo ponto de cada linha d'água = baliza i
        malha_final.append(coluna)
    return malha_final
