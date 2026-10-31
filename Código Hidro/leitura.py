import csv   # Biblioteca padrão do Python para ler arquivos CSV

# ============================================================
# LEITURA DO ARQUIVO DO CASCO
# Este arquivo lê UM único CSV por casco, que traz juntas as balizas, as linhas d'água e as
# cotas, e monta a "tabela_cotas", que é a entrada de construir_linha_baliza() em splines.py.
# Formato esperado do CSV (tabela de cotas):
#   - 1ª linha (cabeçalho): a 1ª célula é só um rótulo (ex.: "x_m"); as demais células são as
#                           alturas z das linhas d'água, em ordem
#   - demais linhas:        a 1ª célula é a posição x da baliza (em metros); as demais são as cotas y
#                           medidas em cada linha d'água (uma coluna por z do cabeçalho)
#   - célula vazia = sem dado naquele cruzamento baliza x linha d'água
#   - linhas em branco e linhas começando com "#" são ignoradas (podem servir de anotação)
# Exemplo:
#   x_m,  0.0,    0.4876, 1.6255
#   0.0,  0.0,    0.0,    0.0
#   3.048,0.0519, 0.2252, 0.3363
# O separador pode ser vírgula ou ponto-e-vírgula, e o decimal pode ser ponto ou vírgula
# (o Excel em português costuma salvar com ponto-e-vírgula e decimal com vírgula).
# ============================================================
#AAAAAAAAA
# ============================================================
# FUNÇÕES AUXILIARES
# Tratam os detalhes do arquivo (número, linhas ignoradas e separador), para a função principal
# trabalhar só com a tabela já limpa.
# ============================================================

def _para_float(texto):
    # Converte o texto de uma célula em número, aceitando decimal com vírgula ("0,4876") ou ponto ("0.4876").
    return float(texto.strip().replace(',', '.'))

def _ler_celulas(caminho_arquivo):
    # Lê o arquivo e devolve a tabela como lista de linhas, cada uma sendo uma lista de células (texto).
    with open(caminho_arquivo, newline='', encoding='utf-8-sig') as arquivo:   # utf-8-sig remove o BOM do Excel, se existir
        linhas = arquivo.read().splitlines()   # todas as linhas do arquivo, sem o "\n" do final

    # descarta linhas em branco (inclusive as que o Excel salva só com separadores, como ";;;;") e comentários
    linhas = [l for l in linhas if l.strip(' ,;\t') != '' and not l.lstrip().startswith('#')]
    if not linhas:
        raise ValueError(f"Arquivo vazio: {caminho_arquivo}")

    # se a primeira linha (cabeçalho) tem ';', esse é o separador (e as vírgulas são decimais); senão, o separador é ','
    separador = ';' if ';' in linhas[0] else ','
    return list(csv.reader(linhas, delimiter=separador))   # separa cada linha em células

# ============================================================
# MONTAGEM DA TABELA DE COTAS
# Para cada baliza, reúne os pares (z, y) medidos e devolve uma lista de tuplas
# (x, lista de z, lista de y), uma por baliza, no formato esperado por
# construir_linha_baliza() em splines.py.
# ============================================================

def ler_tabela_cotas(caminho_arquivo):
    celulas = _ler_celulas(caminho_arquivo)   # tabela do arquivo, já sem linhas em branco/comentários

    # --- CABEÇALHO: alturas z das linhas d'água ---
    colunas_z = {}   # {número da coluna: altura z}; colunas de cabeçalho vazio (sobras do Excel) são ignoradas
    for coluna, texto in enumerate(celulas[0]):
        if coluna == 0 or texto.strip() == '':
            continue                                   # a 1ª célula é só rótulo; células vazias não são linha d'água
        z = _para_float(texto)                         # altura z desta linha d'água
        if z in colunas_z.values():
            raise ValueError(f"Altura z repetida no cabeçalho: {texto.strip()}")   # duas colunas com o mesmo z quebrariam a spline
        colunas_z[coluna] = z

    # --- LEITURA DAS BALIZAS E DAS COTAS ---
    pontos_por_baliza = {}   # {x da baliza: lista de pares (z, y) medidos nessa baliza}
    for numero, linha in enumerate(celulas[1:], start=1):   # "numero" conta só as linhas de dados (1 = primeira baliza)
        try:
            x = _para_float(linha[0])                  # 1ª célula: posição x da baliza, em metros
        except ValueError:
            raise ValueError(f"Linha de dados {numero}: x inválido ('{linha[0]}').")
        if x in pontos_por_baliza:
            raise ValueError(f"Baliza com x = {x} aparece duas vezes.")   # mesma baliza repetida
        pontos_por_baliza[x] = []                      # cria a lista de pontos desta baliza

        for coluna, texto in enumerate(linha):
            if coluna == 0 or texto.strip() == '':
                continue                               # pula o x e as células vazias (sem dado)
            if coluna not in colunas_z:
                raise ValueError(f"Baliza x = {x}: há cota na coluna {coluna + 1}, que não tem altura z no cabeçalho.")
            y = _para_float(texto)                     # cota y (meia-boca) neste cruzamento
            pontos_por_baliza[x].append((colunas_z[coluna], y))   # guarda o par (z, y)

    # --- MONTAGEM DA TABELA FINAL ---
    tabela_cotas = []   # lista de tuplas (x, z_conhecidos, y_conhecidos), uma por baliza
    for x in sorted(pontos_por_baliza.keys()):                          # balizas em ordem crescente de x (0 = proa)
        pontos = sorted(pontos_por_baliza[x])                           # ordena os pares por z crescente (o que splines.py espera)
        if not pontos:
            continue                                                    # baliza sem nenhuma cota: ignora
        z_conhecidos = [ponto[0] for ponto in pontos]                   # separa as alturas z
        y_conhecidos = [ponto[1] for ponto in pontos]                   # separa as cotas y correspondentes
        tabela_cotas.append((x, z_conhecidos, y_conhecidos))
    return tabela_cotas