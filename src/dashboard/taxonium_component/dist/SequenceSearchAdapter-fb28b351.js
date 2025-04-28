import { aV as v, aW as x, aX as F, aY as C, t as n } from "./index-8ee3baa4.js";
import "react";
import "react-dom";
class $ extends v.BaseFeatureDataAdapter {
  async configure() {
    var e;
    const c = await ((e = this.getSubAdapter) === null || e === void 0 ? void 0 : e.call(this, this.getConf("sequenceAdapter")));
    if (!c)
      throw new Error("Error getting subadapter");
    return c.dataAdapter;
  }
  async getRefNames() {
    return (await this.configure()).getRefNames();
  }
  getFeatures(e, c) {
    return x(async (o) => {
      var d;
      const l = await this.configure(), f = 1e3;
      let { start: s, end: r } = e;
      if (s = Math.max(0, s - f), r += f, r < 0 || s > r) {
        o.complete();
        return;
      }
      const u = l.getFeatures({
        ...e,
        start: s,
        end: r
      }, c), m = ((d = (await F.firstValueFrom(u.pipe(C.toArray())))[0]) === null || d === void 0 ? void 0 : d.get("seq")) || "", i = this.getConf("search"), g = this.getConf("searchForward"), w = this.getConf("searchReverse"), A = this.getConf("caseInsensitive"), p = new RegExp(i, "g" + (A ? "i" : ""));
      if (i) {
        if (g) {
          const h = m.matchAll(p);
          for (const a of h) {
            const t = s + (a.index || 0);
            n.doesIntersect2(t, t + i.length, e.start, e.end) && o.next(new n.SimpleFeature({
              uniqueId: `${this.id}-match-${t}-p`,
              refName: e.refName,
              start: t,
              end: t + a[0].length,
              name: a[0],
              strand: 1
            }));
          }
        }
        if (w) {
          const h = n.revcom(m).matchAll(p);
          for (const a of h) {
            const t = r - (a.index || 0);
            n.doesIntersect2(t, t + i.length, e.start, e.end) && o.next(new n.SimpleFeature({
              uniqueId: `${this.id}-match-${t}-n`,
              refName: e.refName,
              start: t - a[0].length,
              name: a[0],
              end: t,
              strand: -1
            }));
          }
        }
      }
      o.complete();
    });
  }
  freeResources() {
  }
}
export {
  $ as default
};
