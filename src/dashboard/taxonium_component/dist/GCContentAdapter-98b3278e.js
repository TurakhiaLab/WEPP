import { aV as A, aW as q, aX as v, aY as x, t as y } from "./index-8ee3baa4.js";
import "react";
import "react-dom";
class F extends A.BaseFeatureDataAdapter {
  constructor() {
    super(...arguments), this.windowSize = 1e3, this.windowDelta = 1e3, this.gcMode = "content";
  }
  async configure() {
    var e;
    const n = await ((e = this.getSubAdapter) === null || e === void 0 ? void 0 : e.call(this, this.getConf("sequenceAdapter")));
    if (!n)
      throw new Error("Error getting subadapter");
    return n.dataAdapter;
  }
  async getRefNames() {
    return (await this.configure()).getRefNames();
  }
  getFeatures(e, n) {
    return q(async (l) => {
      var u;
      const g = await this.configure(), s = this.windowSize === 1 ? 1 : this.windowSize / 2, m = this.windowSize === 1;
      let { start: r, end: o } = e;
      if (r = Math.max(0, r - s), o += s, o < 0 || r > o) {
        l.complete();
        return;
      }
      const S = g.getFeatures({
        ...e,
        start: r,
        end: o
      }, n), p = ((u = (await v.firstValueFrom(S.pipe(x.toArray())))[0]) === null || u === void 0 ? void 0 : u.get("seq")) || "";
      for (let t = s; t < p.length - s; t += this.windowDelta) {
        const i = m ? p[t] : p.slice(t - s, t + s);
        let c = 0, d = 0, h = 0;
        for (let a = 0; a < i.length; a++)
          i[a] === "c" || i[a] === "C" ? c++ : (i[a] === "g" || i[a] === "G") && d++, i[a] !== "N" && h++;
        const f = r;
        let w;
        this.gcMode === "content" ? w = (d + c) / (h || 1) : this.gcMode === "skew" && (w = (d - c) / (d + c || 1)), l.next(new y.SimpleFeature({
          uniqueId: `${this.id}_${f + t}`,
          start: f + t,
          end: f + t + this.windowDelta,
          score: w
        }));
      }
      l.complete();
    });
  }
  /**
   * called to provide a hint that data tied to a certain region
   * will not be needed for the foreseeable future and can be purged
   * from caches, etc
   */
  freeResources() {
  }
}
F.capabilities = ["hasLocalStats"];
export {
  F as default
};
