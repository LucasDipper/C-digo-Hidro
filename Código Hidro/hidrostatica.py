import numpy as np   # Biblioteca numérica: arrays, produto vetorial (np.cross), somas e normas

# ============================================================
# CÁLCULOS HIDROSTÁTICOS A PARTIR DOS PAINÉIS DO CASCO
# Este arquivo recebe os painéis gerados em splines.py e calcula as propriedades hidrostáticas
# (IC "Determinação das propriedades hidrostáticas de unidades flutuantes", cap. 4).
# O arquivo é organizado em camadas:
#   1) PREPARAÇÃO DO CALADO: desloca a origem para a linha d'água e mantém só a parte submersa
#      do casco (corta os painéis que a linha d'água atravessa). Roda UMA vez por calado.
#   2) CARACTERÍSTICAS: uma função para cada propriedade (volume, área molhada, ...). Cada uma
#      recebe a parte submersa já preparada e faz só a sua soma.
#   3) CALADO ÚNICO: prepara um calado e calcula as características dele.
#   4) CURVA: repete o item 3 para uma lista de calados.
# REGRA IMPORTANTE: cada calado é sempre preparado a partir dos painéis ORIGINAIS (nunca a partir
# de um calado já preparado). Por isso a preparação trabalha numa cópia e não altera o original.
# ============================================================

# ============================================================
# ENTRADA: PAINÉIS EM FORMATO DE ARRAY
# splines.py gera os painéis como listas Python. Aqui eles viram um único array numpy, de
# formato (N, 4, 3): N painéis, 4 cantos por painel, 3 coordenadas (x, y, z) por canto.
# Isso permite calcular todos os painéis de uma vez, sem laços, o que deixa o cálculo muito mais rápido.
# ============================================================

def paineis_para_array(paineis):
    # Converte a lista de painéis [[p1, p2, p3, p4], ...] em um array numpy de formato (N, 4, 3).
    # Deve ser chamada UMA vez; o array resultante é o "original" que todo calado usa como base.
    return np.array(paineis, dtype=float)

# ============================================================
# CORTE DE PAINÉIS PELA LINHA D'ÁGUA
# O corte é feito em TRIÂNGULOS: cada painel (quadrilátero) é dividido em dois triângulos, e
# cada triângulo é cortado pelo plano da linha d'água (z = 0, depois da mudança de origem).
# O triângulo é sempre plano, então o corte cobre todos os casos com o mesmo código
# (em vez de tratar um a um os 12 casos da Fig. 4.4 do IC):
#   - 3 vértices submersos: o triângulo fica inteiro
#   - 0 vértices submersos: o triângulo é descartado
#   - 1 ou 2 vértices submersos: o triângulo é cortado e só a parte de baixo é mantida
# A ordem dos vértices (que define o sentido da normal, para fora do casco) é preservada.
# ============================================================

def _dividir_em_triangulos(paineis_array):
    # Divide cada painel (p1, p2, p3, p4) em dois triângulos: (p1, p2, p3) e (p1, p3, p4).
    # A soma dos dois vetores-área é exatamente o vetor-área do painel (metade do vetor de
    # vetor_normal_painel), e o sentido para fora do casco se mantém.
    p1 = paineis_array[:, 0]   # canto 1 de todos os painéis, formato (N, 3)
    p2 = paineis_array[:, 1]   # canto 2 de todos os painéis
    p3 = paineis_array[:, 2]   # canto 3 de todos os painéis
    p4 = paineis_array[:, 3]   # canto 4 de todos os painéis
    triangulos_a = np.stack([p1, p2, p3], axis=1)   # primeiro triângulo de cada painel, formato (N, 3, 3)
    triangulos_b = np.stack([p1, p3, p4], axis=1)   # segundo triângulo de cada painel, formato (N, 3, 3)
    return np.concatenate([triangulos_a, triangulos_b])   # todos os triângulos, formato (2N, 3, 3)

def _cortar_triangulo(triangulo):
    # Corta UM triângulo (3 vértices) pelo plano z = 0 e devolve a parte com z <= 0 (submersa),
    # já dividida em triângulos (uma lista com 0, 1 ou 2 triângulos).
    # Método: percorre as 3 arestas; mantém os vértices submersos e cria um vértice novo
    # em cada aresta que cruza a linha d'água.
    poligono = []   # vértices da parte submersa, em ordem (terá 3 ou 4 vértices)
    for i in range(3):
        atual = triangulo[i]               # vértice atual
        proximo = triangulo[(i + 1) % 3]   # próximo vértice (o último volta para o primeiro)
        atual_submerso = atual[2] <= 0     # True se o vértice atual está na água (ou sobre a linha d'água)
        proximo_submerso = proximo[2] <= 0 # True se o próximo vértice está na água

        if atual_submerso:
            poligono.append(atual)         # vértice submerso fica na parte submersa

        if atual_submerso != proximo_submerso:   # a aresta atravessa a linha d'água
            t = atual[2] / (atual[2] - proximo[2])                  # fração da aresta onde z = 0
            poligono.append(atual + t * (proximo - atual))          # ponto de cruzamento, sobre a linha d'água

    if len(poligono) < 3:
        return []                                                   # nada submerso
    triangulos = [np.array([poligono[0], poligono[1], poligono[2]])]   # primeiro triângulo da parte submersa
    if len(poligono) == 4:
        triangulos.append(np.array([poligono[0], poligono[2], poligono[3]]))   # 4 vértices: sobra um segundo triângulo
    return triangulos

# ============================================================
# 1) PREPARAÇÃO DO CALADO
# Recebe os painéis ORIGINAIS e um calado, e devolve a parte submersa pronta para os cálculos.
# O sistema de coordenadas passa a ser o do IC: origem na linha d'água (z = 0, e z < 0 na água)
# e a meia-nau (x = 0). O eixo x continua crescendo da proa para a popa, como em splines.py;
# por isso o sinal de LCB e LCF deve ser invertido no resultado final, para ficar positivo a vante.
# (Só o z precisa estar na linha d'água para o volume: o painel que fecharia o plano de flutuação
# tem normal vertical e por isso não contribui nos termos x e y; no termo z só contribui zero
# se estiver em z = 0.)
# ============================================================

def preparar_calado(paineis_array, calado, x_meia_nau):
    # paineis_array: painéis ORIGINAIS, formato (N, 4, 3) (saída de paineis_para_array)
    # calado:        altura da linha d'água, medida a partir da quilha (mesmo z da tabela de cotas)
    # x_meia_nau:    posição x da meia-nau, no sistema de splines.py (x = 0 na proa)
    pontos = paineis_array.copy()      # COPIA: o array original nunca é alterado, cada calado parte do original
    pontos[:, :, 0] -= x_meia_nau      # origem em x na meia-nau
    pontos[:, :, 2] -= calado          # origem em z na linha d'água: a água fica em z <= 0

    triangulos = _dividir_em_triangulos(pontos)   # painéis -> triângulos, formato (2N, 3, 3)
    z_max = triangulos[:, :, 2].max(axis=1)       # maior z entre os 3 vértices de cada triângulo
    z_min = triangulos[:, :, 2].min(axis=1)       # menor z entre os 3 vértices de cada triângulo

    inteiros = triangulos[z_max <= 0]                          # triângulos totalmente submersos (ficam como estão)
    parciais = triangulos[(z_max > 0) & (z_min < 0)]           # triângulos atravessados pela linha d'água (precisam de corte)
    # (os que têm z_min >= 0 estão totalmente fora da água e são descartados)

    cortados = []                      # partes submersas dos triângulos atravessados
    for triangulo in parciais:         # só a faixa perto da linha d'água passa por este laço (poucos triângulos)
        cortados.extend(_cortar_triangulo(triangulo))

    if cortados:
        submersos = np.concatenate([inteiros, np.array(cortados)])   # inteiros + partes cortadas
    else:
        submersos = inteiros                                         # nenhum triângulo foi cortado

    a = submersos[:, 0]                # primeiro vértice de cada triângulo submerso
    b = submersos[:, 1]                # segundo vértice
    c = submersos[:, 2]                # terceiro vértice
    vetor_area = 0.5 * np.cross(b - a, c - a)   # vetor-área de cada triângulo: módulo = área, sentido = para fora do casco
    centro = (a + b + c) / 3                    # centro de cada triângulo (média dos 3 vértices; exato para triângulos)

    return {
        'calado': calado,              # calado desta preparação
        'x_meia_nau': x_meia_nau,      # meia-nau usada (no sistema de splines.py)
        'triangulos': submersos,       # vértices dos triângulos submersos, formato (M, 3, 3)
        'vetor_area': vetor_area,      # vetor-área de cada triângulo, formato (M, 3)
        'centro': centro,              # centro de cada triângulo, formato (M, 3)
    }

# ============================================================
# 2) CARACTERÍSTICAS
# Cada função recebe a parte submersa devolvida por preparar_calado() e calcula uma propriedade.
# Aqui entram as mais simples; as demais (plano de flutuação, LCF, TCF, LCB, TCB, KB,
# inércias e BM) serão acrescentadas depois, com a mesma estrutura.
# ============================================================

def termos_do_volume(submerso):
    # Os três termos da Eq. 4.2 do IC: soma de (área vetorial x posição do centro) em cada eixo.
    # Em um casco fechado os três termos são iguais ao volume; compará-los serve de teste de fechamento.
    area = submerso['vetor_area']    # vetor-área de cada triângulo, formato (M, 3)
    centro = submerso['centro']      # centro de cada triângulo, formato (M, 3)
    termo_x = np.sum(area[:, 0] * centro[:, 0])   # soma de A_x * C_x
    termo_y = np.sum(area[:, 1] * centro[:, 1])   # soma de A_y * C_y
    termo_z = np.sum(area[:, 2] * centro[:, 2])   # soma de A_z * C_z
    return termo_x, termo_y, termo_z

def volume(submerso):
    # Volume deslocado (Eq. 4.2): média dos três termos, para reduzir imprecisões numéricas.
    termo_x, termo_y, termo_z = termos_do_volume(submerso)
    return (termo_x + termo_y + termo_z) / 3

def area_molhada(submerso):
    # Área da superfície molhada (Eq. 4.5): soma dos módulos dos vetores-área de todos os triângulos submersos.
    return float(np.sum(np.linalg.norm(submerso['vetor_area'], axis=1)))

def area_plano_flutuacao(submerso):
    # Área do plano de flutuação (Eq. 4.6): menos a soma da componente z dos vetores-área dos triângulos submersos.
    # A componente z é a "sombra" horizontal de cada triângulo; o sinal negativo é porque, no fundo do casco, a normal aponta para baixo.
    area = submerso['vetor_area']        # vetor-área de cada triângulo, formato (M, 3)
    return float(-np.sum(area[:, 2]))    # soma só a coluna z (índice 2) e troca o sinal

# ============================================================
# 3) CALADO ÚNICO
# Prepara UM calado (a partir do original) e calcula todas as características disponíveis.
# Devolve um dicionário, para a curva e o gráfico poderem pegar qualquer característica pelo nome.
# ============================================================

def calcular_hidrostatica(paineis_array, calado, x_meia_nau):
    submerso = preparar_calado(paineis_array, calado, x_meia_nau)   # corte feito UMA vez para este calado
    return {
        'calado': calado,                          # calado deste resultado
        'volume': volume(submerso),                # volume deslocado
        'area_molhada': area_molhada(submerso),    # área da superfície molhada
        'area_plano_flutuacao': area_plano_flutuacao(submerso),     # área do plano de flutuação
    }

# ============================================================
# 4) CURVA
# Repete o cálculo de um calado para uma lista de calados. Calado específico = lista com um só
# calado; curva = lista com vários. Cada calado é preparado a partir do array ORIGINAL
# (paineis_array nunca é alterado), então a ordem dos calados não influencia o resultado.
# ============================================================

def calcular_curva(paineis_array, calados, x_meia_nau):
    # Devolve uma lista de dicionários, um por calado, na mesma ordem de "calados".
    return [calcular_hidrostatica(paineis_array, calado, x_meia_nau) for calado in calados]