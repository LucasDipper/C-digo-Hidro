from splines import fechar_fundo, vetor_normal_painel

malha_teste_fundo = [
    [[0.0, 2.0, 0.0]],   # baliza em x=0, só o ponto da quilha
    [[5.0, 3.0, 0.0]],   # baliza em x=5, só o ponto da quilha
]

paineis = fechar_fundo(malha_teste_fundo)
print("painel gerado:", paineis)

normais = vetor_normal_painel(paineis)
print("componente z da normal:", normais[0].tolist()[2])   # esperado: negativo