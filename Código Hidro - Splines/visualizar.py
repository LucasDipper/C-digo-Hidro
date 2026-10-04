import json
import os
import webbrowser

_TEMPLATE = """<!DOCTYPE html>
<html lang="pt-BR">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>Visualizador da malha do casco</title>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js"></script>
<style>
  :root{ --panel:#1a1f26; --text:#e8ebef; --muted:#9aa3b1; --accent:#5fb4ff; --border:#2a313b; }
  @media (prefers-color-scheme: light){
    :root:not([data-theme="dark"]){ --panel:#ffffff; --text:#1b222b; --muted:#5b6572; --accent:#1f6fd6; --border:#dbe1e8; }
  }
  :root[data-theme="dark"]{ --panel:#1a1f26; --text:#e8ebef; --muted:#9aa3b1; --accent:#5fb4ff; --border:#2a313b; }
  *{box-sizing:border-box;}
  html,body{height:100%;margin:0;background:#11151a;color:var(--text);font-family:ui-sans-serif,system-ui,Segoe UI,Roboto,sans-serif;overflow:hidden;}
  #c{position:fixed;inset:0;display:block;touch-action:none;}
  #panel{position:fixed;top:calc(16px + env(safe-area-inset-top,0px));left:16px;background:var(--panel);
    border:1px solid var(--border);border-radius:10px;padding:14px 16px;max-width:280px;font-size:13px;line-height:1.5;}
  #panel h1{font-size:14px;margin:0 0 8px;font-weight:600;}
  #panel p{margin:0 0 10px;color:var(--muted);}
  #panel input[type=file]{width:100%;font-size:12px;color:var(--text);}
  #panel label{display:flex;gap:8px;align-items:center;margin-top:10px;cursor:pointer;}
  #status{margin-top:8px;color:var(--muted);font-size:12px;}
  #hint{position:fixed;bottom:calc(16px + env(safe-area-inset-bottom,0px));left:16px;color:var(--muted);font-size:12px;}
</style>
</head>
<body>
<canvas id="c"></canvas>
<div id="panel">
  <h1>Malha do casco</h1>
  <p>__TITULO__</p>
  <label><input type="checkbox" id="mirror" checked> Mostrar bombordo (espelhado)</label>
  <label><input type="checkbox" id="paineisToggle" checked> Mostrar superfícies (painéis)</label>
  <p style="margin-top:10px;">Carregar outro arquivo (opcional):</p>
  <input type="file" id="file" accept=".json,application/json">
  <div id="status">Carregando dados embutidos...</div>
</div>
<div id="hint">Arraste para girar · roda do mouse para zoom</div>
<script>
const DADOS_EMBUTIDOS = __DADOS_JSON__;

const canvas = document.getElementById('c');
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x11151a);

const camera = new THREE.PerspectiveCamera(45, innerWidth/innerHeight, 0.1, 100000);
const renderer = new THREE.WebGLRenderer({canvas, antialias:true});
renderer.setPixelRatio(Math.min(devicePixelRatio,2));
renderer.setSize(innerWidth, innerHeight);

let groupMalha = new THREE.Group(); scene.add(groupMalha);
let groupFechamentoLinhas = new THREE.Group(); scene.add(groupFechamentoLinhas);
let groupPaineis = new THREE.Group(); scene.add(groupPaineis);
let groupPerfil = new THREE.Group(); scene.add(groupPerfil);

let target = new THREE.Vector3(0,0,0);
let radius = 200, theta = Math.PI/3, phi = Math.PI/2.5;
let dragging=false, lastX=0, lastY=0;
function updateCamera(){
  const x = target.x + radius*Math.sin(phi)*Math.cos(theta);
  const y = target.y + radius*Math.cos(phi);
  const z = target.z + radius*Math.sin(phi)*Math.sin(theta);
  camera.position.set(x,y,z);
  camera.up.set(0,1,0);
  camera.lookAt(target);
}
canvas.addEventListener('pointerdown', e=>{dragging=true; lastX=e.clientX; lastY=e.clientY;});
addEventListener('pointerup', ()=>dragging=false);
addEventListener('pointermove', e=>{
  if(!dragging) return;
  theta -= (e.clientX-lastX)*0.005;
  phi = Math.min(Math.max(phi-(e.clientY-lastY)*0.005, 0.05), Math.PI-0.05);
  lastX=e.clientX; lastY=e.clientY;
  updateCamera();
});
canvas.addEventListener('wheel', e=>{
  radius = Math.max(1, Math.min(radius*(1+e.deltaY*0.001), 20000));
  updateCamera();
}, {passive:true});
addEventListener('resize', ()=>{
  camera.aspect = innerWidth/innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(innerWidth, innerHeight);
});

const map = p => new THREE.Vector3(p[0], p[2], p[1]);
function corFaixa(z, zmin, zmax){
  const t = zmax>zmin ? (z-zmin)/(zmax-zmin) : 0.5;
  return new THREE.Color().setHSL((1-t)*0.62, 0.65, 0.55);
}
// painéis de fechamento (fundo/pontas) cruzam a linha de centro (2 cantos +y, 2 cantos -y).
// Como cada canto negativo é o espelho exato de um canto positivo no mesmo (x,z), "grudar"
// y em 0 corta o painel exatamente na linha de centro, sem distorcer nada.
function recortarNoCentro(p, mostrarBombordo){
  return mostrarBombordo ? p : [p[0], Math.max(p[1],0), p[2]];
}

function construirMalha(malha, mostrarBombordo){
  scene.remove(groupMalha);
  groupMalha = new THREE.Group();
  scene.add(groupMalha);
  if(!malha) return;

  let zmin=Infinity, zmax=-Infinity;
  malha.flat().forEach(p=>{ if(p){ zmin=Math.min(zmin,p[2]); zmax=Math.max(zmax,p[2]); }});

  const lados = mostrarBombordo ? [1,-1] : [1];
  const pontosGeom = [], coresGeom = [];
  const linhasMat = new THREE.LineBasicMaterial({color:0x5fb4ff, transparent:true, opacity:0.5});
  function flush(pts){ if(pts.length>1) groupMalha.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), linhasMat)); pts.length=0; }

  lados.forEach(lado=>{
    malha.forEach(coluna=>{
      const pts=[];
      coluna.forEach(p=>{
        if(p){
          const v = map([p[0], p[1]*lado, p[2]]);
          pts.push(v);
          pontosGeom.push(v.x,v.y,v.z);
          const c = corFaixa(p[2], zmin, zmax); coresGeom.push(c.r,c.g,c.b);
        } else { flush(pts); }
      });
      flush(pts);
    });
    const numJ = malha[0].length;
    for(let j=0;j<numJ;j++){
      const pts=[];
      malha.forEach(coluna=>{
        const p = coluna[j];
        if(p){ pts.push(map([p[0], p[1]*lado, p[2]])); } else { flush(pts); }
      });
      flush(pts);
    }
  });

  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(pontosGeom,3));
  geo.setAttribute('color', new THREE.Float32BufferAttribute(coresGeom,3));
  groupMalha.add(new THREE.Points(geo, new THREE.PointsMaterial({size:3, vertexColors:true, sizeAttenuation:false})));
}

// desenha o contorno dos painéis de fechamento SEMPRE (independe do toggle de superfícies),
// para nunca deixar um "buraco" visual na costura entre o casco e os fechamentos.
function construirLinhasFechamento(fechamento, mostrarBombordo){
  scene.remove(groupFechamentoLinhas);
  groupFechamentoLinhas = new THREE.Group();
  scene.add(groupFechamentoLinhas);
  if(!fechamento || !fechamento.length) return;
  const mat = new THREE.LineBasicMaterial({color:0x5fb4ff, transparent:true, opacity:0.5});
  fechamento.forEach(painel=>{
    const pts = painel.map(p=>map(recortarNoCentro(p, mostrarBombordo)));
    pts.push(pts[0]);
    groupFechamentoLinhas.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), mat));
  });
}

function adicionarPainel(painel, transform, pos, cores, zmin, zmax){
  const [p1,p2,p3,p4] = painel.map(transform);
  const v1=map(p1), v2=map(p2), v3=map(p3), v4=map(p4);
  const c1=corFaixa(p1[2],zmin,zmax), c2=corFaixa(p2[2],zmin,zmax), c3=corFaixa(p3[2],zmin,zmax), c4=corFaixa(p4[2],zmin,zmax);
  [[v1,v2,v3,c1,c2,c3],[v1,v3,v4,c1,c3,c4]].forEach(([a,b,c,ca,cb,cc])=>{
    pos.push(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z);
    cores.push(ca.r,ca.g,ca.b, cb.r,cb.g,cb.b, cc.r,cc.g,cc.b);
  });
}

function construirPaineis(boreste, bombordo, fechamento, mostrarBombordo, mostrarPaineis){
  scene.remove(groupPaineis);
  groupPaineis = new THREE.Group();
  scene.add(groupPaineis);
  if(!mostrarPaineis) return;

  const todos = [...(boreste||[])];
  if(mostrarBombordo) todos.push(...(bombordo||[]));
  todos.push(...(fechamento||[]));
  if(!todos.length) return;

  let zmin=Infinity, zmax=-Infinity;
  todos.forEach(p=>p.forEach(v=>{ zmin=Math.min(zmin,v[2]); zmax=Math.max(zmax,v[2]); }));

  const pos=[], cores=[];
  (boreste||[]).forEach(p=>adicionarPainel(p, x=>x, pos, cores, zmin, zmax));
  if(mostrarBombordo) (bombordo||[]).forEach(p=>adicionarPainel(p, x=>x, pos, cores, zmin, zmax));
  (fechamento||[]).forEach(p=>adicionarPainel(p, x=>recortarNoCentro(x, mostrarBombordo), pos, cores, zmin, zmax));

  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(pos,3));
  geo.setAttribute('color', new THREE.Float32BufferAttribute(cores,3));
  // polygonOffset evita "z-fighting" (piscar/riscar) entre a superfície e as linhas da malha,
  // que ocupam exatamente as mesmas coordenadas.
  const mat = new THREE.MeshBasicMaterial({vertexColors:true, side:THREE.DoubleSide,
    polygonOffset:true, polygonOffsetFactor:1, polygonOffsetUnits:1});
  groupPaineis.add(new THREE.Mesh(geo, mat));
}

function construirPerfil(perfil) {
  scene.remove(groupPerfil);
  groupPerfil = new THREE.Group();
  scene.add(groupPerfil);
  if(!perfil) return;
  
  // Material da curva (Linha Vermelha)
  const mat = new THREE.LineBasicMaterial({color: 0xff0000, linewidth: 2});

  // silhueta lateral fechada (quilha + proa + convés + popa), vinda do contorno real da malha
  if(perfil.contorno && perfil.contorno.length > 1){
    const pts = perfil.contorno.map(p => map(p));
    groupPerfil.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), mat));
    return;
  }
  
  if(perfil.quilha && perfil.quilha.length > 0){
    const ptsQ = perfil.quilha.map(p => map([p[0], p[1], p[2]]));
    groupPerfil.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(ptsQ), mat));
  }
  
  if(perfil.conves && perfil.conves.length > 0){
    const ptsC = perfil.conves.map(p => map([p[0], p[1], p[2]]));
    groupPerfil.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(ptsC), mat));
  }
  
  // Conecta a proa e a popa (linhas verticais nas extremidades)
  if(perfil.quilha && perfil.conves && perfil.quilha.length > 0 && perfil.conves.length > 0) {
    const p1 = map(perfil.quilha[0]); 
    const c1 = map(perfil.conves[0]);
    groupPerfil.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([p1, c1]), mat));
    
    const p2 = map(perfil.quilha[perfil.quilha.length-1]); 
    const c2 = map(perfil.conves[perfil.conves.length-1]);
    groupPerfil.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([p2, c2]), mat));
  }
}

function ajustarCamera(){
  const box = new THREE.Box3().setFromObject(scene);
  if(box.isEmpty()) return;
  const size = box.getSize(new THREE.Vector3());
  target.copy(box.getCenter(new THREE.Vector3()));
  radius = Math.max(size.x,size.y,size.z) * 1.6 || 100;
  updateCamera();
}

function carregar(){
  const mostrarBombordo = document.getElementById('mirror').checked;
  const mostrarPaineis = document.getElementById('paineisToggle').checked;
  construirMalha(window.ultimaMalha, mostrarBombordo);
  construirLinhasFechamento(window.ultimoFechamento, mostrarBombordo);
  construirPaineis(window.ultimoBoreste, window.ultimoBombordo, window.ultimoFechamento, mostrarBombordo, mostrarPaineis);
  
  construirPerfil(window.ultimoPerfil);
  
  ajustarCamera();
}

function aplicarDados(dados, rotulo){
  const malha = Array.isArray(dados) ? dados : dados.malha_final;
  window.ultimaMalha = malha;
  window.ultimoBoreste = Array.isArray(dados) ? [] : (dados.boreste || []);
  window.ultimoBombordo = Array.isArray(dados) ? [] : (dados.bombordo || []);
  window.ultimoFechamento = Array.isArray(dados) ? [] : (dados.fechamento || []);
  
  window.ultimoPerfil = Array.isArray(dados) ? null : dados.perfil_diametral;
  
  carregar();
  const totalPontos = malha ? malha.flat().filter(Boolean).length : 0;
  const totalPaineis = window.ultimoBoreste.length + window.ultimoBombordo.length + window.ultimoFechamento.length;
  document.getElementById('status').textContent =
    `${malha ? malha.length : 0} balizas, ${totalPontos} pontos, ${totalPaineis} painéis (${rotulo}).`;
}

document.getElementById('file').addEventListener('change', e=>{
  const f = e.target.files[0];
  if(!f) return;
  const reader = new FileReader();
  reader.onload = ()=>{
    try{ aplicarDados(JSON.parse(reader.result), 'arquivo carregado'); }
    catch(err){ document.getElementById('status').textContent = 'Arquivo inválido: ' + err.message; }
  };
  reader.readAsText(f);
});
document.getElementById('mirror').addEventListener('change', carregar);
document.getElementById('paineisToggle').addEventListener('change', carregar);

updateCamera();
aplicarDados(DADOS_EMBUTIDOS, 'carregados automaticamente');
(function animate(){ requestAnimationFrame(animate); renderer.render(scene, camera); })();
</script>
</body>
</html>
"""


# ============================================================
# GERAÇÃO DO VISUALIZADOR 3D (arquivo HTML)
# Esta função recebe a malha, os painéis e o perfil lateral calculados pelo pipeline, embute tudo
# em um modelo de página HTML (_TEMPLATE) e salva o resultado em um arquivo que abre no navegador.
# Novidade desta versão: o argumento perfil_diametral, com a silhueta lateral do casco.
# Depende de, no início do arquivo visualizar.py (não incluídos aqui): os módulos
# json, os e webbrowser, e a variável _TEMPLATE, que guarda o HTML/JavaScript do visualizador
# com os marcadores __DADOS_JSON__ e __TITULO__.
# ============================================================

def gerar_visualizacao_3d(malha_final, paineis_boreste=None, paineis_bombordo=None, paineis_fechamento=None, perfil_diametral=None,
                           caminho_saida="visualizacao_casco.html",
                           titulo="Gerado automaticamente pelo pipeline Python.", abrir_navegador=True):
    """
    Gera um arquivo HTML autocontido com o visualizador 3D do casco, com os
    dados já embutidos (sem precisar de upload manual). Os painéis são
    passados separados por categoria:
      - paineis_boreste: painéis do casco original (y positivo)
      - paineis_bombordo: painéis espelhados (y negativo)
      - paineis_fechamento: painéis do fundo e das pontas, que cruzam a
        linha de centro (2 cantos de cada lado)
    Essa separação permite ligar/desligar o bombordo corretamente, inclusive
    recortando os painéis de fechamento exatamente na linha de centro.
    """
    # --- MONTAGEM DOS DADOS ---
    # Reúne tudo o que o visualizador precisa em um único dicionário.
    # "x or []" (ou "x or {}") troca None por uma lista (ou dicionário) vazia, para o visualizador
    # sempre receber um valor do tipo esperado.
    dados = {
        "malha_final": malha_final,                # malha de pontos [x, y, z]
        "boreste": paineis_boreste or [],          # painéis do lado de boreste
        "bombordo": paineis_bombordo or [],        # painéis do lado de bombordo (espelhados)
        "fechamento": paineis_fechamento or [],    # painéis de fundo, proa e popa
        "perfil_diametral": perfil_diametral or {}   # silhueta lateral (plano de simetria); no main: {"contorno": contorno_perfil(malha_final)}
    }

    # --- PREENCHIMENTO DO MODELO HTML ---
    # json.dumps converte o dicionário em texto JSON, que é inserido no lugar do marcador
    # __DADOS_JSON__ do modelo; assim os dados ficam dentro do próprio HTML.
    html = _TEMPLATE.replace("__DADOS_JSON__", json.dumps(dados))
    html = html.replace("__TITULO__", titulo)      # coloca o título no marcador correspondente

    # --- GRAVAÇÃO DO ARQUIVO ---
    caminho_absoluto = os.path.abspath(caminho_saida)   # transforma o caminho em absoluto (necessário para abrir no navegador)
    with open(caminho_absoluto, "w", encoding="utf-8") as f:   # "w" cria o arquivo (ou sobrescreve se já existir)
        f.write(html)

    # --- ABERTURA NO NAVEGADOR ---
    if abrir_navegador:
        webbrowser.open(f"file://{caminho_absoluto}")   # abre o HTML recém-gerado no navegador padrão

    print(f"Visualização salva em: {caminho_absoluto}")   # informa onde o arquivo foi salvo
    return caminho_absoluto                                # devolve o caminho, caso o chamador queira usá-lo