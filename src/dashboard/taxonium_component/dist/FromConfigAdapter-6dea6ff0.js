import { aV as c, l as i, ar as m, aW as l } from "./index-8ee3baa4.js";
import "react";
import "react-dom";
class o extends c.BaseFeatureDataAdapter {
  constructor(e, s, t) {
    super(e, s, t);
    const a = i.readConfObject(e, "features");
    this.features = o.makeFeatures(a || []);
  }
  static makeFeatures(e) {
    const s = /* @__PURE__ */ new Map();
    for (let t = 0; t < e.length; t += 1)
      if (e[t]) {
        const a = this.makeFeature(e[t]), r = a.get("refName");
        let u = s.get(r);
        u || (u = [], s.set(r, u)), u.push(a);
      }
    for (const t of s.values())
      t.sort((a, r) => a.get("start") - r.get("start"));
    return s;
  }
  static makeFeature(e) {
    return new m(e);
  }
  async getRefNames() {
    return [...this.features.keys()];
  }
  async getRefNameAliases() {
    return [...this.features.values()].map((e) => ({
      refName: e[0].get("refName"),
      aliases: e[0].get("aliases")
    }));
  }
  getFeatures(e) {
    const { refName: s, start: t, end: a } = e;
    return l(async (r) => {
      const u = this.features.get(s) || [];
      for (let f = 0; f < u.length; f++) {
        const n = u[f];
        n.get("end") > t && n.get("start") < a && r.next(n);
      }
      r.complete();
    });
  }
  freeResources() {
  }
}
export {
  o as default
};
