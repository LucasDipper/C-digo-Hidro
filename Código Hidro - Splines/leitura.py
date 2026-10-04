import csv   # módulo padrão do Python para ler arquivos CSV - não precisa instalar nada

# ============================================================
# LEITURA DE ARQUIVOS CSV
# Este arquivo lê os três CSVs de cada embarcação (balizas, linhas d'água e cotas) e
# monta a "tabela_cotas", que é a entrada de construir_linha_baliza() em splines.py.
# Colunas esperadas em cada CSV:
#   - balizas:       "indice" e "x_m"                    (posição x de cada baliza, em metros)
#   - linhas d'água: "indice" e "z"                      (altura z de cada linha d'água)
#   - cotas:         "baliza", "linha_dagua" e "cota_y"  (y medido em cada cruzamento baliza x linha d'água)
# ============================================================

# ============================================================
# FUNÇÃO AUXILIAR: cria o leitor de CSV
# Usada pelas três funções de leitura, para tratar do mesmo jeito arquivos com
# separador vírgula ou ponto-e-vírgula.
# ============================================================

def _abrir_leitor(caminho_arquivo, arquivo_aberto):
    # Recebe um arquivo JÁ aberto e devolve um leitor (DictReader) que entrega cada linha como dicionário
    # O encoding utf-8-sig, que remove o BOM (marca invisível no início do arquivo, comum em CSVs
    # salvos pelo Excel) se ele existir e não faz nada se não existir, é definido no open()
    # de quem chama esta função.
    # csv.Sniffer detecta sozinho se o separador é vírgula ou ponto-e-vírgula.
    amostra = arquivo_aberto.read(2048)   # lê só o começo do arquivo para o Sniffer analisar
    arquivo_aberto.seek(0)                # volta ao início, para o DictReader ler o arquivo todo
    try:
        dialeto = csv.Sniffer().sniff(amostra, delimiters=',;')   # tenta descobrir o separador (',' ou ';')
    except csv.Error:
        dialeto = csv.excel   # padrão: vírgula (usado quando o Sniffer não consegue decidir)
    return csv.DictReader(arquivo_aberto, dialect=dialeto)   # a 1ª linha do CSV vira o nome das colunas


# ============================================================
# LEITURA DAS BALIZAS E DAS LINHAS D'ÁGUA
# Cada uma devolve um dicionário {índice: valor}, para depois ligar os índices da tabela de cotas
# às posições reais (x das balizas, z das linhas d'água).
# ============================================================

def ler_balizas(caminho_arquivo):
    # Lê balizas.csv e devolve um dicionário {indice: x_m},
    # para depois conseguir traduzir "baliza número tal" em "posição x tal".
    balizas_x = {}   # dicionário de resultado
    # o "with open(...) as arquivo" garante que o arquivo é fechado sozinho
    # no final, mesmo se der erro no meio - não precisa chamar arquivo.close()
    with open(caminho_arquivo, newline='', encoding='utf-8-sig') as arquivo:
        # csv.DictReader lê cada linha como um dicionário, usando o cabeçalho
        # (a primeira linha do CSV) como as chaves - assim você acessa por
        # nome de coluna (linha['x_m']) em vez de por posição (linha[2])
        leitor = _abrir_leitor(caminho_arquivo, arquivo)
        for linha in leitor:
            indice = int(linha['indice'])   # tudo que vem do CSV chega como texto - precisa converter
            x = float(linha['x_m'])         # posição longitudinal da baliza, em metros
            balizas_x[indice] = x
    return balizas_x

def ler_linhas_dagua(caminho_arquivo):
    # Mesma lógica de ler_balizas, mas para linhas_dagua.csv -> {indice: z}
    linhas_z = {}    # dicionário de resultado
    with open(caminho_arquivo, newline='', encoding='utf-8-sig') as arquivo:
        leitor = _abrir_leitor(caminho_arquivo, arquivo)
        for linha in leitor:
            indice = int(linha['indice'])   # número da linha d'água (inteiro)
            z = float(linha['z'])           # altura da linha d'água
            linhas_z[indice] = z
    return linhas_z

# ============================================================
# MONTAGEM DA TABELA DE COTAS
# Combina os três arquivos. Para cada baliza, reúne os pares (z, y) medidos e devolve
# uma lista de tuplas (x, lista de z, lista de y), uma por baliza, no formato esperado por
# construir_linha_baliza() em splines.py.
# ============================================================

def ler_tabela_cotas(caminho_balizas, caminho_linhas_dagua, caminho_cotas):
    # Junta os três arquivos e monta a tabela no formato que o resto do
    # código já entende: uma lista de tuplas (x_baliza, [z conhecidos], [y conhecidos]).
    balizas_x = ler_balizas(caminho_balizas)             # {índice da baliza: x}
    linhas_z = ler_linhas_dagua(caminho_linhas_dagua)    # {índice da linha d'água: z}

    # agrupa os pontos (z, y) de cada baliza num dicionário temporário:
    # {indice_baliza: [(z1, y1), (z2, y2), ...]}
    pontos_por_baliza = {}

    # --- LEITURA DAS COTAS ---
    with open(caminho_cotas, newline='', encoding='utf-8-sig') as arquivo:
        leitor = _abrir_leitor(caminho_cotas, arquivo)
        for linha in leitor:
            indice_baliza = int(linha['baliza'])                # a qual baliza pertence esta cota
            indice_linha_dagua = int(linha['linha_dagua'])      # em qual linha d'água foi medida
            y = float(linha['cota_y'])                          # valor da cota y (meia-boca)

            z = linhas_z[indice_linha_dagua]   # traduz o índice da linha d'água para a altura real em metros

            if indice_baliza not in pontos_por_baliza:
                pontos_por_baliza[indice_baliza] = []   # primeira vez que essa baliza aparece: cria a lista
            pontos_por_baliza[indice_baliza].append((z, y))   # guarda o par (z, y)

    # --- MONTAGEM DA TABELA FINAL ---
    # monta a lista final, uma tupla por baliza, com os pontos ordenados por z crescente
    # (necessário porque cotas.csv nem sempre lista as linhas d'água em ordem, e a spline exige isso)
    tabela_cotas = []
    for indice_baliza in sorted(pontos_por_baliza.keys()):   # balizas em ordem crescente de índice
        x = balizas_x[indice_baliza]                          # posição x da baliza
        pontos = sorted(pontos_por_baliza[indice_baliza])   # ordena as tuplas (z,y) pelo z (primeiro item)
        z_conhecidos = [ponto[0] for ponto in pontos]        # separa as alturas z
        y_conhecidos = [ponto[1] for ponto in pontos]        # separa as cotas y correspondentes
        tabela_cotas.append((x, z_conhecidos, y_conhecidos))

    return tabela_cotas