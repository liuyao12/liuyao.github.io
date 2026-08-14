  function displayedNodePoint(kind,index){
    const point=pointsFor(kind)[index], members=groupIndices(kind,index), center=toScreen(point,state.zView,zPlot);
    if(members.length<=1) return center;
    const k=members.indexOf(index), r=Math.min(7,3.4+members.length*.55), a=2*Math.PI*k/members.length-Math.PI/2;
    return {x:center.x+r*Math.cos(a),y:center.y+r*Math.sin(a)};
  }

  function mergeGroups(kind,index,targetIndex){
    const groups=groupsFor(kind), points=pointsFor(kind), a=groups[index], b=groups[targetIndex];
    if(a===b)return;
    const id=Math.min(a,b), target=cloneC(points[targetIndex]);
    for(let i=0;i<groups.length;i++) if(groups[i]===a||groups[i]===b){groups[i]=id;points[i]=cloneC(target);}
    if(kind==='z')setFromRoots(points,true);else setFromW({keepC:true});
  }

  function unsnapGroup(kind,index){
    const members=groupIndices(kind,index); if(members.length<=1)return false;
    const points=pointsFor(kind), groups=groupsFor(kind), center=cloneC(points[index]);
    const pixelUnit=2*state.zView.half/Math.min(zPlot.clientWidth,zPlot.clientHeight), r=pixelUnit*(members.length===2?10:12);
    let next=freshGroupId(groups);
    members.forEach((i,k)=>{points[i]=add(center,scale(expi(2*Math.PI*k/members.length-Math.PI/2),r));groups[i]=next++;});
    if(kind==='z')setFromRoots(points,true);else setFromW({keepC:true});
    return true;
  }

  function nearestSnapTarget(kind,index){
    const points=pointsFor(kind), groups=groupsFor(kind), own=groups[index], p=toScreen(points[index],state.zView,zPlot);
    let best=-1,bestD=18;
    for(let i=0;i<points.length;i++)if(groups[i]!==own){const q=toScreen(points[i],state.zView,zPlot),d=Math.hypot(p.x-q.x,p.y-q.y);if(d<bestD){bestD=d;best=i;}}
    return best;
  }

  // ---------- SVG render ----------
  function renderZ(full=true){
    const zd=dims(zPlot);zSvg.setAttribute('viewBox',`0 0 ${zd.w} ${zd.h}`);clearSvg(zSvg);drawGrid(zSvg,state.zView,zPlot);drawLandscape(full);
    const showDisk=document.getElementById('diskToggle').checked;
    state.flowData=computeFlowData(full);

    const hull=convexHull(state.Z);
    if(hull.length>=2){let d='';hull.forEach((z,k)=>{const p=toScreen(z,state.zView,zPlot);d+=(k?'L':'M')+p.x.toFixed(2)+','+p.y.toFixed(2);});if(hull.length>=3)d+='Z';zSvg.appendChild(svgEl('path',{d,class:'focus-boundary'}));}

    if(showDisk){const disk=state.flowData.disk,p=toScreen(disk.center,state.zView,zPlot),pixelsPerUnit=Math.min(zd.w,zd.h)/(2*state.zView.half);zSvg.appendChild(svgEl('circle',{cx:p.x,cy:p.y,r:disk.radius*pixelsPerUnit,class:'enclosing-disk'}));zSvg.appendChild(svgEl('line',{x1:p.x-5,y1:p.y,x2:p.x+5,y2:p.y,class:'disk-center'}));zSvg.appendChild(svgEl('line',{x1:p.x,y1:p.y-5,x2:p.x,y2:p.y+5,class:'disk-center'}));}

    for(const branch of state.flowData.ridges){if(branch.points.length<2)continue;let d='';branch.points.forEach((z,k)=>{const p=toScreen(z,state.zView,zPlot);d+=(k?'L':'M')+p.x.toFixed(2)+','+p.y.toFixed(2);});zSvg.appendChild(svgEl('path',{d,class:'flow ridge'}));}
    const selectedIds=new Set(state.flowData.best.filter(x=>x.branch).map(x=>x.branch.id));
    for(const branch of state.flowData.branches){if(branch.points.length<2)continue;let d='';branch.points.forEach((z,k)=>{const p=toScreen(z,state.zView,zPlot);d+=(k?'L':'M')+p.x.toFixed(2)+','+p.y.toFixed(2);});const cls=branch.reachedRoot?(selectedIds.has(branch.id)?'flow down-selected':'flow down-secondary'):'flow unresolved';zSvg.appendChild(svgEl('path',{d,class:cls}));}

    for(const kind of ['z','w']){
      const points=pointsFor(kind),groups=groupsFor(kind),seen=new Set();
      points.forEach((point,i)=>{
        const gid=groups[i],members=groupIndices(kind,i),center=toScreen(point,state.zView,zPlot);
        if(members.length>1&&!seen.has(gid)){seen.add(gid);zSvg.appendChild(svgEl('circle',{cx:center.x,cy:center.y,r:10+Math.min(5,members.length),class:'stack-ring'}));}
        const p=displayedNodePoint(kind,i),node=svgEl('circle',{cx:p.x,cy:p.y,r:7.1,class:`node ${kind==='z'?'root-node':'critical-node'}`});
        node.addEventListener('pointerdown',e=>startDrag(e,kind,i));zSvg.appendChild(node);
      });
    }
    drawCompactValueControl(zSvg);
  }


  function renderFlowBars(){
    const data=state.flowData || computeFlowData();
    const bars=document.getElementById('flowBars');
    const axis=document.getElementById('flowAxis');
    const summary=document.getElementById('flowSummary');
    const ordered=data.best.slice().sort((a,b)=>{
      if(!Number.isFinite(a.ratio) && !Number.isFinite(b.ratio)) return a.rootIndex-b.rootIndex;
      if(!Number.isFinite(a.ratio)) return -1;
      if(!Number.isFinite(b.ratio)) return 1;
      return b.ratio-a.ratio;
    });
    const finite=ordered.filter(x=>Number.isFinite(x.ratio));
    const maxRatio=finite.length?Math.max(...finite.map(x=>x.ratio)):1;
    const scaleMax=Math.max(1.08,Math.ceil(maxRatio*10+1)/10);
    const thresholdPct=Math.min(100,100/scaleMax);

    bars.innerHTML=ordered.map(item=>{
      const ok=Number.isFinite(item.ratio);
      const pct=ok?Math.min(100,100*item.ratio/scaleMax):0;
      const over=ok&&item.ratio>1+1e-6;
      const value=ok?item.ratio.toFixed(3):'—';
      const title=ok
        ? `arc length ${item.length.toPrecision(5)}; normalized ${item.ratio.toPrecision(5)}`
        : 'No descending saddle-to-root branch was resolved numerically.';
      return `<div class="flow-bar-row" title="${title}">
        <span class="flow-bar-root"></span>
        <span class="flow-bar-track">
          <span class="flow-bar-fill${over?' over':''}" style="width:${pct.toFixed(3)}%"></span>
          <span class="flow-bar-threshold" style="left:${thresholdPct.toFixed(3)}%"></span>
        </span>
        <span class="flow-bar-value${over?' over':''}${ok?'':' unresolved'}">${value}</span>
      </div>`;
    }).join('');

    const ticks=[0,1];
    if(scaleMax>1.12) ticks.push(scaleMax);
    axis.innerHTML=ticks.map(v=>`<span style="left:${(100*v/scaleMax).toFixed(3)}%">${v===scaleMax&&v!==1?v.toFixed(1):v}</span>`).join('');

    const unresolved=ordered.filter(x=>!Number.isFinite(x.ratio)).length;
    const worst=finite.length?Math.max(...finite.map(x=>x.ratio)):NaN;
    const verdict=Number.isFinite(worst)
      ? (worst<=1+1e-6 ? `<b>All traced lengths ≤ R.</b> Largest: ${worst.toFixed(4)}.` : `<b>${finite.filter(x=>x.ratio>1+1e-6).length} traced length${finite.filter(x=>x.ratio>1+1e-6).length===1?'':'s'} exceed R.</b> Largest: ${worst.toFixed(4)}.`)
      : 'No paths resolved.';
    summary.innerHTML=verdict+(unresolved?` ${unresolved} root${unresolved===1?'':'s'} unresolved numerically.`:'');
  }

  function renderText(){
    renderFlowBars();
    const data=state.flowData||computeFlowData(false), finite=data.best.filter(x=>Number.isFinite(x.ratio));
    const worst=finite.length?Math.max(...finite.map(x=>x.ratio)):NaN;
    document.getElementById('compactStatus').textContent=`R ${data.R.toFixed(3)}${Number.isFinite(worst)?` · max ${worst.toFixed(3)}`:''}`;
  }
  function escapeHtml(s){return s.replace(/[&<>]/g,ch=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[ch]));}

  function scheduleRender(full=false){
    state.fullRender=state.fullRender||full;
    if(state.renderQueued) return;
    state.renderQueued=true;
    requestAnimationFrame(()=>{
      state.renderQueued=false; const f=state.fullRender; state.fullRender=false;
      renderZ(f);renderText();
    });
  }

  // ---------- dragging ----------
  function startDrag(e,kind,index){
    e.preventDefault();e.stopPropagation();
    const p=evtPoint(e,zPlot);
    state.dragging={kind,index,pointerId:e.pointerId,startX:p.x,startY:p.y,moved:false};
    if(kind==='z'||kind==='w'){
      const members=groupIndices(kind,index), points=pointsFor(kind);
      state.dragging.members=members;
      state.dragging.initial=members.map(i=>cloneC(points[i]));
      state.dragging.startWorld=fromScreen(p.x,p.y,state.zView,zPlot);
    }else if(kind==='c'){
      state.dragging.valueScale=compactValueScale();
      state.dragging.widget=compactWidgetGeometry();
    }
    e.target.setPointerCapture?.(e.pointerId);
  }

  window.addEventListener('pointermove',e=>{
    const d=state.dragging;if(!d||e.pointerId!==d.pointerId)return;
    const p=evtPoint(e,zPlot);if(Math.hypot(p.x-d.startX,p.y-d.startY)>3)d.moved=true;
    if(d.kind==='z'||d.kind==='w'){
      const target=fromScreen(p.x,p.y,state.zView,zPlot),delta=sub(target,d.startWorld),points=pointsFor(d.kind);
      d.members.forEach((idx,k)=>points[idx]=add(d.initial[k],delta));
      if(d.kind==='z')setFromRoots(points,true);else setFromW({keepC:true});
    }else if(d.kind==='c'){
      setC(compactFromScreen(p.x,p.y,d.valueScale,d.widget));
    }
    scheduleRender(false);
  });

  window.addEventListener('pointerup',e=>{
    const d=state.dragging;if(!d||e.pointerId!==d.pointerId)return;
    if(d.kind==='z'||d.kind==='w'){
      if(!d.moved)unsnapGroup(d.kind,d.index);
      else{const target=nearestSnapTarget(d.kind,d.index);if(target>=0)mergeGroups(d.kind,d.index,target);}
    }
    state.dragging=null;scheduleRender(true);
  });
  window.addEventListener('pointercancel',()=>{state.dragging=null;scheduleRender(true);});

  function wheelZoom(e,view,el){
    e.preventDefault();
    const p=evtPoint(e,el), before=fromScreen(p.x,p.y,view,el);
    const factor=Math.exp(Math.sign(e.deltaY)*.14); view.half=Math.max(.15,Math.min(100,view.half*factor));
    const after=fromScreen(p.x,p.y,view,el); view.cx+=before.re-after.re; view.cy+=before.im-after.im;
    scheduleRender(true);
  }
  zPlot.addEventListener('wheel',e=>wheelZoom(e,state.zView,zPlot),{passive:false});

  // ---------- controls ----------
  document.getElementById('degreeSelect').addEventListener('change',e=>reset(+e.target.value));
  document.getElementById('resetBtn').addEventListener('click',()=>reset(state.n));
  document.getElementById('randomBtn').addEventListener('click',randomize);
  document.getElementById('fitBtn').addEventListener('click',fitViews);
  for(const id of ['landscapeToggle','diskToggle'])document.getElementById(id).addEventListener('change',()=>scheduleRender(true));
  new ResizeObserver(()=>scheduleRender(true)).observe(zPlot);

  reset(5);
