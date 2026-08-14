  // ---------- landscape and flow rendering ----------
  function resizeCanvas(canvas,el){
    const dpr=1, w=Math.max(1,Math.floor(el.clientWidth*dpr)), h=Math.max(1,Math.floor(el.clientHeight*dpr));
    if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;}
    return {w,h,dpr};
  }

  // Subdued domain coloring uses arg(P).  A separate red/blue tint, computed
  // from the local critical expansion of P, marks uphill and downhill sectors
  // without inserting an artificial hue wheel around the critical point.
  function hsvToRgb(h,s,v){
    h=((h%1)+1)%1; s=Math.max(0,Math.min(1,s)); v=Math.max(0,Math.min(1,v));
    const q=h*6, i=Math.floor(q), f=q-i;
    const p=v*(1-s), a=v*(1-s*f), b=v*(1-s*(1-f));
    const rgb=[[v,b,p],[a,v,p],[p,v,b],[p,a,v],[b,p,v],[v,p,a]][i%6];
    return rgb.map(x=>Math.round(255*x));
  }

  function periodicLine(x,width){
    const d=Math.abs(x-Math.round(x));
    return Math.exp(-0.5*(d/width)**2);
  }

  function mixRgb(a,b,t){
    t=Math.max(0,Math.min(1,t));
    return a.map((x,i)=>Math.round(x+(b[i]-x)*t));
  }

  function criticalLocalExpansion(w){
    let q=polyDerivative(state.dP), order=2;
    const scaleRef=Math.max(1,abs(polyEval(state.P,w)));
    while(order<=state.n){
      const value=polyEval(q,w);
      if(abs(value)>1e-8*scaleRef) return {order,coefficient:scale(value,1/factorial(order))};
      q=polyDerivative(q);order++;
    }
    return null;
  }
  function factorial(n){let r=1;for(let k=2;k<=n;k++)r*=k;return r;}

  function prepareSaddleData(){
    const disk=smallestEnclosingCircle(state.Z);
    const R=Math.max(disk.radius,1e-6);
    return state.W.map((w,index)=>{
      const pw=polyEval(state.P,w), local=criticalLocalExpansion(w);
      if(abs(pw)<1e-11 || !local) return null;
      let nearestRoot=Infinity, nearestSaddle=Infinity;
      for(const z of state.Z) nearestRoot=Math.min(nearestRoot,dist(w,z));
      for(let k=0;k<state.W.length;k++) if(k!==index) nearestSaddle=Math.min(nearestSaddle,dist(w,state.W[k]));
      const radius=Math.max(.28*R,Math.min(.78*R,.72*nearestRoot,.55*nearestSaddle));
      return {w,pw,A:div(local.coefficient,pw),order:local.order,radius,logLevel:Math.log(Math.max(1e-300,abs(pw)))};
    }).filter(Boolean);
  }

  function prepareFocusData(){
    const feather=Math.max(.025,state.zView.half*.10),hull=convexHull(state.Z);
    if(hull.length<3)return {mode:'disk',feather,disk:smallestEnclosingCircle(state.Z)};
    return {mode:'hull',feather,hull};
  }


  function applyFocusDimming(rgb,z,focus){
    if(!focus || focus.mode==='none') return rgb;
    const signed=focus.mode==='disk'
      ? focus.disk.radius-dist(z,focus.disk.center)
      : signedDistanceToHull(z,focus.hull);
    const outside=smoothstep(0,focus.feather,-signed);
    if(outside<=0) return rgb;
    const luminance=Math.round(.2126*rgb[0]+.7152*rgb[1]+.0722*rgb[2]);
    let muted=mixRgb(rgb,[luminance,luminance,luminance],.72*outside);
    muted=mixRgb(muted,[7,16,24],.58*outside);
    return muted;
  }

  function domainColor(z,value,derivative,saddles,focus){
    const m=Math.max(1e-300,abs(value));
    const phase=(arg(value)+Math.PI)/(2*Math.PI);
    const log2m=Math.log2(m);

    // Ordinary domain coloring remains visible, but deliberately subdued.
    const major=periodicLine(log2m,.055), minor=periodicLine(log2m-.5,.08);
    let v=.46+.12*Math.tanh(log2m/4.5)+.18*major+.045*minor;
    v=Math.max(.14,Math.min(.88,v));
    let rgb=hsvToRgb(phase,.46,v);

    // The nearest saddle tints the actual landscape: red uphill, blue downhill.
    // Close to the saddle we use its first nonzero local term, so all alternating
    // sectors remain visible, including at a multiple critical point.
    if(saddles.length){
      let best=null,bestScaled=Infinity;
      for(const q of saddles){
        const delta=sub(z,q.w), scaled=abs(delta)/Math.max(q.radius,1e-9);
        if(scaled<bestScaled){bestScaled=scaled;best={q,delta};}
      }
      if(best){
        const {q,delta}=best;
        let power=C(1,0);for(let k=0;k<q.order;k++)power=mul(power,delta);
        const leading=mul(q.A,power);
        const denom=abs(q.A)*Math.max(Math.pow(abs(delta),q.order),1e-20);
        const local=Math.max(-1,Math.min(1,leading.re/denom));
        const exact=Math.tanh(5.2*(Math.log(m)-q.logLevel));
        const near=Math.exp(-1.8*bestScaled*bestScaled);
        const signed=near*local+(1-near)*exact;
        const up=[236,73,78], down=[55,104,222];
        const target=signed>=0?up:down;
        const strength=(.28+.62/(1+Math.pow(bestScaled/1.18,3)))*(.42+.58*Math.abs(signed));
        rgb=mixRgb(rgb,target,Math.min(.90,strength));
        const ridge=Math.exp(-.5*Math.pow(signed/.075,2))*Math.exp(-.7*bestScaled*bestScaled);
        rgb=mixRgb(rgb,[240,244,246],.22*ridge);
      }
    }
    return applyFocusDimming(rgb,z,focus);
  }

  function drawLandscape(highQuality=true){
    const show=document.getElementById('landscapeToggle').checked;
    const {w,h}=resizeCanvas(landscapeCanvas,zPlot);
    const ctx=landscapeCanvas.getContext('2d');
    ctx.clearRect(0,0,w,h);
    if(!show) return;

    // Redraw at a coarser resolution while a point is moving, then sharpen
    // on pointer release.  This keeps the color field synchronized with P.
    const side=Math.min(w,h);
    const N=highQuality
      ? Math.max(130,Math.min(230,Math.round(side/3.3)))
      : Math.max(70,Math.min(105,Math.round(side/7.0)));
    const off=document.createElement('canvas'); off.width=N; off.height=N;
    const octx=off.getContext('2d',{alpha:false});
    const img=octx.createImageData(N,N);
    const saddles=prepareSaddleData();
    const focus=prepareFocusData();

    for(let j=0;j<N;j++) for(let i=0;i<N;i++){
      const z=fromScreen((i+.5)/N*zPlot.clientWidth,(j+.5)/N*zPlot.clientHeight,state.zView,zPlot);
      const value=polyEval(state.P,z);
      const derivative=polyEval(state.dP,z);
      const [r,g,b]=domainColor(z,value,derivative,saddles,focus);
      const k=4*(j*N+i);
      img.data[k]=r; img.data[k+1]=g; img.data[k+2]=b; img.data[k+3]=255;
    }
    octx.putImageData(img,0,0);
    ctx.imageSmoothingEnabled=true;
    ctx.drawImage(off,0,0,w,h);
  }


  function traceSeparatrices(w, criticalIndex, disk, highQuality=true){
    const pw=polyEval(state.P,w), local=criticalLocalExpansion(w);
    if(abs(pw)<1e-9 || !local) return [];
    const B=mul(conj(pw),local.coefficient), order=local.order;
    const geomScale=Math.max(disk.radius,1e-4);
    const eps=Math.max(geomScale*.0014,1e-7);
    const step=Math.max(geomScale*(highQuality?.0046:.012),2e-7);
    const hitRadius=Math.max(step*2.6,geomScale*3e-4);
    const escapeRadius=Math.max(8*geomScale,4*state.zView.half);
    const maxSteps=highQuality?1900:520;
    const paths=[];

    for(const type of ['down','up']){
      const fieldSign=type==='down'?-1:1;
      for(let arm=0;arm<order;arm++){
        const theta=((type==='down'?Math.PI:0)-arg(B)+2*Math.PI*arm)/order;
        let z=add(w,scale(expi(theta),eps));
        const pts=[cloneC(w),cloneC(z)];
        let length=eps,reachedRoot=false,rootIndex=null,reachedCritical=false;
        let lastStored=cloneC(z);
        for(let k=0;k<maxSteps;k++){
          if(type==='down'){
            let nearest=Infinity,nearestIndex=-1;
            for(let i=0;i<state.Z.length;i++){const d=dist(z,state.Z[i]);if(d<nearest){nearest=d;nearestIndex=i;}}
            if(nearest<hitRadius){length+=nearest;z=cloneC(state.Z[nearestIndex]);if(dist(z,lastStored)>1e-12)pts.push(cloneC(z));reachedRoot=true;rootIndex=nearestIndex;break;}
          }else{
            for(let j=0;j<state.W.length;j++) if(j!==criticalIndex && dist(z,state.W[j])<hitRadius*1.4){z=cloneC(state.W[j]);pts.push(cloneC(z));reachedCritical=true;break;}
            if(reachedCritical)break;
          }
          const p=polyEval(state.P,z),dp=polyEval(state.dP,z);
          let v=scale(mul(p,conj(dp)),fieldSign);
          const vm=abs(v); if(vm<1e-15||!finiteC(v))break; v=scale(v,1/vm);
          const mid=add(z,scale(v,step*.5));
          let v2=scale(mul(polyEval(state.P,mid),conj(polyEval(state.dP,mid))),fieldSign);
          const v2m=abs(v2);if(v2m<1e-15||!finiteC(v2))break;v2=scale(v2,1/v2m);
          const next=add(z,scale(v2,step));length+=dist(z,next);z=next;
          if(k%(highQuality?2:1)===0){pts.push(cloneC(z));lastStored=cloneC(z);}
          if(dist(z,w)>escapeRadius||!finiteC(z))break;
        }
        paths.push({id:`${criticalIndex}:${type}:${arm}`,type,criticalIndex,rootIndex,reachedRoot,reachedCritical,length,points:pts});
      }
    }
    return paths;
  }

  function computeFlowData(highQuality=true){
    const disk=smallestEnclosingCircle(state.Z),branches=[],ridges=[];
    const best=state.Z.map((z,rootIndex)=>({rootIndex,length:Infinity,branch:null,direct:false}));
    const coincidenceTol=Math.max(1e-8,disk.radius*1e-7);
    for(let i=0;i<state.Z.length;i++) for(let j=0;j<state.W.length;j++) if(dist(state.Z[i],state.W[j])<=coincidenceTol){best[i]={rootIndex:i,length:0,branch:null,direct:true,criticalIndex:j};break;}
    state.W.forEach((w,j)=>{
      for(const branch of traceSeparatrices(w,j,disk,highQuality)){
        if(branch.type==='up'){ridges.push(branch);continue;}
        branches.push(branch);
        if(branch.reachedRoot&&branch.length<best[branch.rootIndex].length)best[branch.rootIndex]={rootIndex:branch.rootIndex,length:branch.length,branch,direct:false,criticalIndex:j};
      }
    });
    const R=disk.radius;
    for(const item of best)item.ratio=Number.isFinite(item.length)?(R>1e-14?item.length/R:0):Infinity;
    return {disk,branches,ridges,best,R};
  }

  function drawCompactValueControl(svg){
    const g=compactWidgetGeometry(), values=criticalValues(), s=compactValueScale(values);
    const group=svgEl('g',{'aria-label':'compactified value-plane control'});
    const bg=svgEl('circle',{cx:g.cx,cy:g.cy,r:g.radius,class:'compact-bg'});
    group.appendChild(bg);
    for(const frac of [.50,.75]) group.appendChild(svgEl('circle',{cx:g.cx,cy:g.cy,r:g.inner*frac,class:'compact-ring'}));
    group.appendChild(svgEl('line',{x1:g.cx-g.inner,y1:g.cy,x2:g.cx+g.inner,y2:g.cy,class:'compact-axis'}));
    group.appendChild(svgEl('line',{x1:g.cx,y1:g.cy-g.inner,x2:g.cx,y2:g.cy+g.inner,class:'compact-axis'}));
    const caption=svgEl('text',{x:g.cx,y:g.cy+g.radius-8,'text-anchor':'middle',class:'compact-caption'});
    caption.textContent='c (compactified)'; group.appendChild(caption);
    const inf=svgEl('text',{x:g.cx,y:g.cy-g.radius+11,'text-anchor':'middle',class:'compact-caption'});
    inf.textContent='∞'; group.appendChild(inf);

    for(const v of values){
      const p=compactToScreen(v,s,g);
      group.appendChild(svgEl('circle',{cx:p.x,cy:p.y,r:3.4,class:'critical-value-node'}));
    }
    const cp=compactToScreen(state.c,s,g);
    const node=svgEl('circle',{cx:cp.x,cy:cp.y,r:7.2,class:'node parameter-node'});
    node.addEventListener('pointerdown',e=>startDrag(e,'c',0));
    group.appendChild(node);

    group.addEventListener('pointerdown',e=>{
      startDrag(e,'c',0);
      const p=evtPoint(e,zPlot);
      setC(compactFromScreen(p.x,p.y,s,g));
      scheduleRender(false);
    });
    svg.appendChild(group);
  }

