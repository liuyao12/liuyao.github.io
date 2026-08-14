  // ---------- state ----------
  const state={
    n:5,
    Z:[], W:[], P:[], dP:[], F:[], c:C(0,0),
    zGroups:[], wGroups:[],
    dragging:null,
    zView:{cx:0,cy:0,half:2.5},
    flowData:null,
    renderQueued:false,
    fullRender:true
  };

  const singletonGroups=n=>Array.from({length:n},(_,i)=>i);
  function groupsFor(kind){ return kind==='z'?state.zGroups:state.wGroups; }
  function pointsFor(kind){ return kind==='z'?state.Z:state.W; }
  function setGroupsFor(kind,g){ if(kind==='z')state.zGroups=g;else state.wGroups=g; }
  function groupIndices(kind,index){ const g=groupsFor(kind), id=g[index]; return g.map((x,i)=>x===id?i:-1).filter(i=>i>=0); }
  function freshGroupId(g){ return g.length?Math.max(...g)+1:0; }
  function resetGroups(kind){ setGroupsFor(kind,singletonGroups(pointsFor(kind).length)); }

  const defaults={
    3:[C(-1.35,-.45),C(.05,1.25),C(1.25,-.55)],
    4:[C(-1.45,-.55),C(-.55,1.15),C(.72,.92),C(1.35,-.72)],
    5:[C(-1.55,-.58),C(-.72,1.12),C(.38,1.30),C(1.52,.15),C(.36,-1.30)],
    6:[C(-1.62,-.44),C(-1.08,.91),C(-.18,1.42),C(.92,1.05),C(1.55,-.18),C(.23,-1.38)],
    7:[C(-1.65,-.25),C(-1.25,.85),C(-.44,1.45),C(.48,1.34),C(1.33,.72),C(1.57,-.51),C(.08,-1.48)],
    8:[C(-1.7,-.2),C(-1.35,.82),C(-.65,1.42),C(.24,1.53),C(1.08,1.02),C(1.62,.25),C(1.22,-1.02),C(-.18,-1.48)]
  };

  function setFromRoots(roots, preserveZGroups=false){
    const oldW=state.W;
    state.Z=roots.map(cloneC);
    state.n=state.Z.length;
    if(!preserveZGroups || state.zGroups.length!==state.n) state.zGroups=singletonGroups(state.n);
    state.P=polyFromRoots(state.Z);
    state.dP=polyDerivative(state.P);
    const rawW=solvePolynomial(state.dP, oldW);
    state.W=matchPoints(oldW,rawW);
    state.wGroups=singletonGroups(state.W.length);
    state.F=polyIntegralZero(state.dP);
    state.c=polyEval(state.F,state.Z[0]);
  }

  function setFromW({keepC=true}={}){
    const oldZ=state.Z.map(cloneC);
    state.dP=polyScale(polyFromRoots(state.W),state.n);
    state.F=polyIntegralZero(state.dP);
    if(!keepC) state.c=C(0,0);
    state.P=state.F.map(cloneC);
    state.P[0]=sub(state.P[0],state.c);
    const rawZ=solvePolynomial(state.P,oldZ);
    state.Z=matchPoints(oldZ,rawZ);
    state.zGroups=singletonGroups(state.Z.length);
    if(state.wGroups.length!==state.W.length) state.wGroups=singletonGroups(state.W.length);
  }

  function setC(c){
    const oldZ=state.Z.map(cloneC);
    state.c=cloneC(c);
    state.P=state.F.map(cloneC);
    state.P[0]=sub(state.P[0],state.c);
    const rawZ=solvePolynomial(state.P,oldZ);
    state.Z=matchPoints(oldZ,rawZ);
    state.zGroups=singletonGroups(state.Z.length);
  }

  function reset(n=state.n){
    state.n=n;
    state.zGroups=[]; state.wGroups=[];
    setFromRoots((defaults[n]||defaults[5]).map(cloneC));
    fitViews();
  }


  function randomize(){
    const roots=[];
    for(let k=0;k<state.n;k++){
      const t=2*Math.PI*k/state.n + (Math.random()-.5)*.45;
      const r=.75+Math.random()*.9;
      roots.push(C(r*Math.cos(t)+(Math.random()-.5)*.25,r*Math.sin(t)+(Math.random()-.5)*.25));
    }
    setFromRoots(roots);
    fitViews();
  }

  // ---------- plot transforms ----------
  const zPlot=document.getElementById('zPlot');
  const zSvg=document.getElementById('zSvg');
  const landscapeCanvas=document.getElementById('landscapeCanvas');
  const NS='http://www.w3.org/2000/svg';
  function dims(el){ return {w:el.clientWidth,h:el.clientHeight}; }
  function toScreen(z,view,el){ const {w,h}=dims(el); const s=Math.min(w,h)/(2*view.half); return {x:w/2+(z.re-view.cx)*s,y:h/2-(z.im-view.cy)*s}; }
  function fromScreen(x,y,view,el){ const {w,h}=dims(el); const s=Math.min(w,h)/(2*view.half); return C(view.cx+(x-w/2)/s,view.cy-(y-h/2)/s); }
  function evtPoint(e,el){ const r=el.getBoundingClientRect(); return {x:e.clientX-r.left,y:e.clientY-r.top}; }
  function svgEl(name,attrs={}){ const el=document.createElementNS(NS,name); for(const [k,v] of Object.entries(attrs)) el.setAttribute(k,v); return el; }
  function clearSvg(svg){ while(svg.firstChild) svg.removeChild(svg.firstChild); }

  function fitViewToPoints(view,pts,minHalf=1){
    if(!pts.length) return;
    let minx=Infinity,maxx=-Infinity,miny=Infinity,maxy=-Infinity;
    for(const z of pts){minx=Math.min(minx,z.re);maxx=Math.max(maxx,z.re);miny=Math.min(miny,z.im);maxy=Math.max(maxy,z.im);}
    view.cx=(minx+maxx)/2; view.cy=(miny+maxy)/2;
    view.half=Math.max(minHalf,(maxx-minx)/2,(maxy-miny)/2)*1.38;
  }
  function fitViews(){
    const disk=smallestEnclosingCircle(state.Z);
    const diskExtrema=[
      C(disk.center.re+disk.radius,disk.center.im),C(disk.center.re-disk.radius,disk.center.im),
      C(disk.center.re,disk.center.im+disk.radius),C(disk.center.re,disk.center.im-disk.radius)
    ];
    fitViewToPoints(state.zView,[...state.Z,...state.W,...diskExtrema],1);
    scheduleRender(true);
  }

  function criticalValues(){ return state.W.map(w=>polyEval(state.F,w)); }

  function compactValueScale(values=criticalValues()){
    const radii=values.map(abs).filter(r=>r>1e-10).sort((a,b)=>a-b);
    if(!radii.length){
      const disk=smallestEnclosingCircle(state.Z);
      return Math.max(1e-4,Math.pow(Math.max(.35,disk.radius),state.n));
    }
    const median=radii[Math.floor((radii.length-1)/2)], max=radii[radii.length-1];
    return Math.max(1e-6,.62*median+.38*max);
  }

  function compactWidgetGeometry(){
    const {w,h}=dims(zPlot);
    const radius=Math.max(58,Math.min(82,.125*Math.min(w,h)));
    return {cx:w-radius-14,cy:radius+14,radius,inner:radius-10};
  }

  function compactMap(v,s){
    const r=abs(v);
    return r<1e-30?C(0,0):scale(v,1/(s+r));
  }

  function compactToScreen(v,s,g){
    const u=compactMap(v,s);
    return {x:g.cx+g.inner*u.re,y:g.cy-g.inner*u.im};
  }

  function compactFromScreen(x,y,s,g){
    let u=C((x-g.cx)/g.inner,-(y-g.cy)/g.inner), r=abs(u);
    const cap=.978;
    if(r>cap){u=scale(u,cap/r);r=cap;}
    return r<1e-12?C(0,0):scale(u,s/(1-r));
  }

  function drawGrid(svg,view,el){
    const {w,h}=dims(el);
    const step=niceStep(2*view.half/5);
    const xmin=view.cx-view.half*w/Math.min(w,h), xmax=view.cx+view.half*w/Math.min(w,h);
    const ymin=view.cy-view.half*h/Math.min(w,h), ymax=view.cy+view.half*h/Math.min(w,h);
    for(let x=Math.ceil(xmin/step)*step;x<=xmax;x+=step){
      const p=toScreen(C(x,0),view,el);
      svg.appendChild(svgEl('line',{x1:p.x,y1:0,x2:p.x,y2:h,class:x===0?'axis-line':'grid-line'}));
    }
    for(let y=Math.ceil(ymin/step)*step;y<=ymax;y+=step){
      const p=toScreen(C(0,y),view,el);
      svg.appendChild(svgEl('line',{x1:0,y1:p.y,x2:w,y2:p.y,class:y===0?'axis-line':'grid-line'}));
    }
  }
  function niceStep(x){ const p=10**Math.floor(Math.log10(x)); const m=x/p; return (m<1.5?1:m<3.5?2:m<7.5?5:10)*p; }

