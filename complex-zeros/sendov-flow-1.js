"use strict";

  // ---------- complex arithmetic ----------
  const C = (re=0, im=0) => ({re, im});
  const add = (a,b) => C(a.re+b.re, a.im+b.im);
  const sub = (a,b) => C(a.re-b.re, a.im-b.im);
  const mul = (a,b) => C(a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re);
  const div = (a,b) => { const d=b.re*b.re+b.im*b.im || 1e-30; return C((a.re*b.re+a.im*b.im)/d,(a.im*b.re-a.re*b.im)/d); };
  const scale = (a,s) => C(a.re*s,a.im*s);
  const neg = a => C(-a.re,-a.im);
  const conj = a => C(a.re,-a.im);
  const abs2 = a => a.re*a.re+a.im*a.im;
  const abs = a => Math.hypot(a.re,a.im);
  const arg = a => Math.atan2(a.im,a.re);
  const expi = t => C(Math.cos(t),Math.sin(t));
  const cloneC = a => C(a.re,a.im);
  const dist = (a,b) => abs(sub(a,b));
  const finiteC = a => Number.isFinite(a.re) && Number.isFinite(a.im);

  // coefficients are ascending: a[0] + a[1] z + ...
  const polyEval = (a,z) => { let y=C(0,0); for(let k=a.length-1;k>=0;k--) y=add(mul(y,z),a[k]); return y; };
  const polyDerivative = a => a.length<=1 ? [C(0,0)] : a.slice(1).map((v,k)=>scale(v,k+1));
  const polyIntegralZero = a => [C(0,0), ...a.map((v,k)=>scale(v,1/(k+1)))];
  const polyScale = (a,s) => a.map(v=>scale(v,s));
  const polyFromRoots = roots => {
    let a=[C(1,0)];
    for(const r of roots){
      const b=Array(a.length+1).fill(0).map(()=>C(0,0));
      for(let k=0;k<a.length;k++){
        b[k]=add(b[k],neg(mul(a[k],r)));
        b[k+1]=add(b[k+1],a[k]);
      }
      a=b;
    }
    return a;
  };

  function solvePolynomial(coeffs, seeds=null){
    let a=coeffs.map(cloneC);
    while(a.length>1 && abs(a[a.length-1])<1e-14) a.pop();
    const n=a.length-1;
    if(n<=0) return [];
    const lead=a[n];
    a=a.map(x=>div(x,lead));
    const da=polyDerivative(a);
    let roots;
    if(seeds && seeds.length===n && seeds.every(finiteC)){
      roots=seeds.map((z,k)=>add(z,scale(expi(1.618*k+.3),1e-5)));
    } else {
      let R=1;
      for(let k=0;k<n;k++) R=Math.max(R,1+abs(a[k]));
      roots=Array.from({length:n},(_,k)=>scale(expi(2*Math.PI*(k+.37)/n),R));
    }
    for(let iter=0;iter<180;iter++){
      let maxCorr=0;
      const next=roots.map((z,i)=>{
        const p=polyEval(a,z), dp=polyEval(da,z);
        let newton=abs(dp)>1e-16 ? div(p,dp) : scale(expi(iter+i),1e-4);
        let s=C(0,0);
        for(let j=0;j<n;j++) if(j!==i){
          let d=sub(z,roots[j]);
          if(abs(d)<1e-12) d=add(d,scale(expi(i+j+iter),1e-8));
          s=add(s,div(C(1,0),d));
        }
        let den=sub(C(1,0),mul(newton,s));
        let corr=abs(den)>1e-12 ? div(newton,den) : newton;
        if(!finiteC(corr) || abs(corr)>10) corr=scale(expi(i+iter*.2),1e-3);
        maxCorr=Math.max(maxCorr,abs(corr));
        return sub(z,corr);
      });
      roots=next;
      if(maxCorr<1e-12) break;
    }
    return roots;
  }

  function matchPoints(oldPts,newPts){
    const n=newPts.length;
    if(!oldPts || oldPts.length!==n) return newPts;
    // exact minimum assignment by subset DP; n <= 8.
    const size=1<<n;
    const dp=Array(size).fill(Infinity), parent=Array(size).fill(null);
    dp[0]=0;
    for(let mask=0;mask<size;mask++){
      const i=popcount(mask);
      if(i>=n || !Number.isFinite(dp[mask])) continue;
      for(let j=0;j<n;j++) if(!(mask&(1<<j))){
        const m2=mask|(1<<j), val=dp[mask]+abs2(sub(oldPts[i],newPts[j]));
        if(val<dp[m2]){dp[m2]=val; parent[m2]=[mask,j];}
      }
    }
    const assign=Array(n); let mask=size-1;
    for(let i=n-1;i>=0;i--){ const [pm,j]=parent[mask]; assign[i]=newPts[j]; mask=pm; }
    return assign;
  }
  function popcount(x){ let c=0; while(x){x&=x-1;c++;} return c; }

  // The Euclidean smallest enclosing circle.  For a finite nonempty set it is
  // unique and is supported by either two antipodal boundary points or three
  // boundary points.  Since n <= 8, exhaustive candidate checking is simple
  // and numerically more transparent than a randomized Welzl implementation.
  function smallestEnclosingCircle(points){
    if(!points.length) return {center:C(0,0), radius:0, support:[]};
    const candidates=[];
    for(let i=0;i<points.length;i++) candidates.push({center:cloneC(points[i]),radius:0,support:[i]});
    for(let i=0;i<points.length;i++) for(let j=i+1;j<points.length;j++){
      const center=scale(add(points[i],points[j]),.5);
      candidates.push({center,radius:dist(center,points[i]),support:[i,j]});
    }
    for(let i=0;i<points.length;i++) for(let j=i+1;j<points.length;j++) for(let k=j+1;k<points.length;k++){
      const a=points[i],b=points[j],c=points[k];
      const d=2*(a.re*(b.im-c.im)+b.re*(c.im-a.im)+c.re*(a.im-b.im));
      if(Math.abs(d)<1e-12) continue;
      const aa=abs2(a),bb=abs2(b),cc=abs2(c);
      const center=C(
        (aa*(b.im-c.im)+bb*(c.im-a.im)+cc*(a.im-b.im))/d,
        (aa*(c.re-b.re)+bb*(a.re-c.re)+cc*(b.re-a.re))/d
      );
      candidates.push({center,radius:dist(center,a),support:[i,j,k]});
    }
    let best=null;
    for(const q of candidates){
      const tol=1e-9*Math.max(1,q.radius);
      if(points.every(z=>dist(z,q.center)<=q.radius+tol) && (!best || q.radius<best.radius-1e-10)) best=q;
    }
    // The pair candidates guarantee that a valid circle is found.
    return best || {center:cloneC(points[0]),radius:0,support:[0]};
  }

  function cross(o,a,b){
    return (a.re-o.re)*(b.im-o.im)-(a.im-o.im)*(b.re-o.re);
  }

  function convexHull(points){
    if(points.length<=1) return points.map(cloneC);
    const pts=points.map(cloneC).sort((a,b)=>a.re-b.re || a.im-b.im);
    const unique=[];
    for(const q of pts) if(!unique.length || dist(q,unique[unique.length-1])>1e-12) unique.push(q);
    if(unique.length<=2) return unique;
    const lower=[];
    for(const q of unique){
      while(lower.length>=2 && cross(lower[lower.length-2],lower[lower.length-1],q)<=1e-12) lower.pop();
      lower.push(q);
    }
    const upper=[];
    for(let i=unique.length-1;i>=0;i--){
      const q=unique[i];
      while(upper.length>=2 && cross(upper[upper.length-2],upper[upper.length-1],q)<=1e-12) upper.pop();
      upper.push(q);
    }
    lower.pop(); upper.pop();
    return lower.concat(upper);
  }

  function pointSegmentDistance(p,a,b){
    const ab=sub(b,a), ap=sub(p,a), den=abs2(ab);
    if(den<1e-24) return dist(p,a);
    const t=Math.max(0,Math.min(1,(ap.re*ab.re+ap.im*ab.im)/den));
    return dist(p,add(a,scale(ab,t)));
  }

  function signedDistanceToHull(p,hull){
    if(hull.length===0) return -Infinity;
    if(hull.length===1) return -dist(p,hull[0]);
    if(hull.length===2) return -pointSegmentDistance(p,hull[0],hull[1]);
    let inside=true, d=Infinity;
    for(let i=0;i<hull.length;i++){
      const a=hull[i],b=hull[(i+1)%hull.length];
      if(cross(a,b,p)<-1e-11) inside=false;
      d=Math.min(d,pointSegmentDistance(p,a,b));
    }
    return inside?d:-d;
  }

  function smoothstep(edge0,edge1,x){
    if(edge0===edge1) return x<edge0?0:1;
    let t=(x-edge0)/(edge1-edge0); t=Math.max(0,Math.min(1,t));
    return t*t*(3-2*t);
  }

  function formatC(z, digits=3){
    const r=Math.abs(z.re)<.5*10**(-digits)?0:z.re;
    const i=Math.abs(z.im)<.5*10**(-digits)?0:z.im;
    if(i===0) return r.toFixed(digits);
    if(r===0) return `${i.toFixed(digits)}i`;
    return `${r.toFixed(digits)} ${i>=0?'+':'−'} ${Math.abs(i).toFixed(digits)}i`;
  }
  function polynomialString(a){
    const parts=[];
    for(let k=a.length-1;k>=0;k--){
      if(abs(a[k])<1e-7) continue;
      const coef=formatC(a[k],3);
      const term=k===0?`(${coef})`:k===1?`(${coef})z`:`(${coef})z^${k}`;
      parts.push(term);
    }
    return parts.join(' + ').replaceAll('+ (−','− (');
  }

