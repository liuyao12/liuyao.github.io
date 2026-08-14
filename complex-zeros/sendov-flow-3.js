  // ---------- contour and flow rendering ----------
  function resizeCanvas(canvas,el){
    const dpr=1, w=Math.max(1,Math.floor(el.clientWidth*dpr)), h=Math.max(1,Math.floor(el.clientHeight*dpr));
    if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;}
    return {w,h,dpr};
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

  function quantile(sorted,q){
    if(!sorted.length)return 0;
    const x=(sorted.length-1)*q, i=Math.floor(x), t=x-i;
    return sorted[i]*(1-t)+sorted[Math.min(sorted.length-1,i+1)]*t;
  }
  function niceContourStep(span,target=11){
    const raw=Math.max(span/target,1e-6), p=10**Math.floor(Math.log10(raw)), m=raw/p;
    return (m<1.5?1:m<3?2:m<7?5:10)*p;
  }
  function autoContourLevels(values){
    const a=values.filter(Number.isFinite).slice().sort((x,y)=>x-y);
    if(a.length<2)return [];
    let lo=quantile(a,.08), hi=quantile(a,.92);
    if(!(hi>lo)){lo=a[0];hi=a[a.length-1];}
    const step=niceContourStep(hi-lo,12), first=Math.ceil(lo/step)*step, out=[];
    for(let x=first;x<=hi+step*1e-7 && out.length<22;x+=step)out.push(x);
    return out;
  }
  function dedupeLevels(levels,tol=.015){
    const a=levels.filter(Number.isFinite).sort((x,y)=>x-y), out=[];
    for(const x of a)if(!out.length||Math.abs(x-out[out.length-1])>tol)out.push(x);
    return out;
  }
  function edgePoint(edge,x,y,a,b,c,d,level){
    const interp=(v0,v1)=>{
      const den=v1-v0;
      return Math.max(0,Math.min(1,Math.abs(den)<1e-14?.5:(level-v0)/den));
    };
    if(edge===0){const t=interp(a,b);return [x+t,y];}
    if(edge===1){const t=interp(b,c);return [x+1,y+t];}
    if(edge===2){const t=interp(d,c);return [x+t,y+1];}
    const t=interp(a,d);return [x,y+t];
  }
  function marchLevel(ctx,grid,N,level,sx,sy){
    ctx.beginPath();
    for(let y=0;y<N;y++)for(let x=0;x<N;x++){
      const row=N+1, a=grid[y*row+x], b=grid[y*row+x+1], c=grid[(y+1)*row+x+1], d=grid[(y+1)*row+x];
      if(!Number.isFinite(a+b+c+d))continue;
      const mask=(a>level?1:0)|(b>level?2:0)|(c>level?4:0)|(d>level?8:0);
      if(mask===0||mask===15)continue;
      const crossings=[];
      if((a>level)!==(b>level))crossings.push(0);
      if((b>level)!==(c>level))crossings.push(1);
      if((d>level)!==(c>level))crossings.push(2);
      if((a>level)!==(d>level))crossings.push(3);
      const segment=(e0,e1)=>{
        const p=edgePoint(e0,x,y,a,b,c,d,level),q=edgePoint(e1,x,y,a,b,c,d,level);
        ctx.moveTo(p[0]*sx,p[1]*sy);ctx.lineTo(q[0]*sx,q[1]*sy);
      };
      if(crossings.length===2){segment(crossings[0],crossings[1]);continue;}
      if(crossings.length!==4)continue;
      const center=(a+b+c+d)/4, centerHigh=center>level;
      if(mask===5){
        if(centerHigh){segment(0,1);segment(2,3);}else{segment(0,3);segment(1,2);}
      }else if(mask===10){
        if(centerHigh){segment(0,3);segment(1,2);}else{segment(0,1);segment(2,3);}
      }else{
        segment(crossings[0],crossings[1]);segment(crossings[2],crossings[3]);
      }
    }
    ctx.stroke();
  }

  function drawLandscape(highQuality=true){
    const show=document.getElementById('landscapeToggle').checked;
    const {w,h}=resizeCanvas(landscapeCanvas,zPlot), ctx=landscapeCanvas.getContext('2d');
    ctx.clearRect(0,0,w,h);
    if(!show)return;

    // Contours of log|P|.  The special levels log|P(w)| pass through the
    // critical points and reveal the saddle crossings directly.
    const side=Math.min(w,h);
    const N=highQuality
      ? Math.max(105,Math.min(185,Math.round(side/4.1)))
      : Math.max(55,Math.min(82,Math.round(side/8.3)));
    const row=N+1, grid=new Float64Array(row*row), samples=[];
    for(let j=0;j<=N;j++)for(let i=0;i<=N;i++){
      const z=fromScreen(i/N*zPlot.clientWidth,j/N*zPlot.clientHeight,state.zView,zPlot);
      const m=Math.max(1e-300,abs(polyEval(state.P,z)));
      const v=Math.max(-40,Math.min(40,Math.log(m)));
      grid[j*row+i]=v;samples.push(v);
    }
    const criticalLevels=dedupeLevels(state.W.map(w=>{
      const m=abs(polyEval(state.P,w));
      return m>1e-12?Math.log(m):NaN;
    }),highQuality?.012:.025);
    const regular=autoContourLevels(samples).filter(x=>criticalLevels.every(c=>Math.abs(x-c)>.045));
    const sx=w/N, sy=h/N;

    ctx.save();
    ctx.lineJoin='round';ctx.lineCap='round';
    ctx.strokeStyle='rgba(197,216,226,.19)';ctx.lineWidth=highQuality?.85:.75;
    for(const level of regular)marchLevel(ctx,grid,N,level,sx,sy);
    ctx.strokeStyle='rgba(239,246,249,.72)';ctx.lineWidth=highQuality?1.55:1.25;
    for(const level of criticalLevels)marchLevel(ctx,grid,N,level,sx,sy);
    ctx.restore();
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
