import { t as T, aV as D, aZ as y, b9 as x, aW as j } from "./index-8ee3baa4.js";
import { I as $ } from "./main.esm-53240cd4.js";
import "react";
import "react-dom";
function B(c) {
  return c[0] === 31 && c[1] === 139 && c[2] === 8;
}
function w(c, e, s, r) {
  const t = c.split("	"), o = t[s ? 3 : 0], h = +t[s ? 4 : 1], n = +t[s ? 5 : 2], u = +t[s ? 0 : 3], i = +t[s ? 1 : 4], a = +t[s ? 2 : 5], d = t[6], l = +t[7], v = b(t[8]), f = b(t[9]), p = t.slice(9), m = r ? Object.fromEntries(r.slice(9).map((g, F) => [g, p[F]])) : p;
  return new T.SimpleFeature({
    start: h,
    end: n,
    refName: o,
    strand: v,
    name: d,
    ...m,
    score: l,
    uniqueId: e,
    mate: { refName: u, start: i, end: a, strand: f }
  });
}
function b(c) {
  return c === "+" ? 1 : c === "-" ? -1 : c === "." ? 0 : void 0;
}
class E extends D.BaseFeatureDataAdapter {
  constructor() {
    super(...arguments), this.intervalTrees = {};
  }
  async loadDataP(e = {}) {
    const s = this.pluginManager, r = this.getConf("bedpeLocation"), t = await y.openLocation(r, s).readFile(e), o = B(t) ? await x(t) : t;
    if (o.length > 536870888)
      throw new Error("Data exceeds maximum string length (512MB)");
    const n = new TextDecoder("utf8", { fatal: !0 }).decode(o).split(/\n|\r\n|\r/).filter((f) => !!f), u = [];
    let i = 0;
    for (; i < n.length && n[i].startsWith("#"); i++)
      u.push(n[i]);
    const a = u.join(`
`), d = {}, l = {};
    for (; i < n.length; i++) {
      const f = n[i], p = f.split("	"), m = p[0], g = p[3];
      d[m] || (d[m] = []), l[g] || (l[g] = []), d[m].push(f), l[g].push(f);
    }
    const v = this.getConf("columnNames");
    return {
      header: a,
      feats1: d,
      feats2: l,
      columnNames: v
    };
  }
  async loadData(e = {}) {
    return this.bedpeFeatures || (this.bedpeFeatures = this.loadDataP(e).catch((s) => {
      throw this.bedpeFeatures = void 0, s;
    })), this.bedpeFeatures;
  }
  async getRefNames(e = {}) {
    const { feats1: s, feats2: r } = await this.loadData(e);
    return [.../* @__PURE__ */ new Set([...Object.keys(s), ...Object.keys(r)])];
  }
  async getHeader(e = {}) {
    const { header: s } = await this.loadData(e);
    return s;
  }
  async getNames() {
    const { header: e, columnNames: s } = await this.loadData();
    if (s.length)
      return s;
    const r = e.split(/\n|\r\n|\r/).filter((o) => !!o), t = r[r.length - 1];
    return t != null && t.includes("	") ? t.slice(1).split("	").map((o) => o.trim()) : void 0;
  }
  async loadFeatureTreeP(e) {
    const { feats1: s, feats2: r } = await this.loadData(), t = s[e], o = r[e], h = await this.getNames(), n = new $(), u = t == null ? void 0 : t.map((a, d) => {
      const l = `${this.id}-${e}-${d}`;
      return w(a, l, !1, h);
    }), i = o == null ? void 0 : o.map((a, d) => {
      const l = `${this.id}-${e}-${d}`;
      return w(a, l, !0, h);
    });
    for (const a of u)
      n.insert([a.get("start"), a.get("end")], a);
    for (const a of i)
      n.insert([a.get("start"), a.get("end")], a);
    return n;
  }
  async loadFeatureTree(e) {
    return this.intervalTrees[e] || (this.intervalTrees[e] = this.loadFeatureTreeP(e).catch((s) => {
      throw this.intervalTrees[e] = void 0, s;
    })), this.intervalTrees[e];
  }
  getFeatures(e, s = {}) {
    return j(async (r) => {
      const { start: t, end: o, refName: h } = e, n = await this.loadFeatureTree(h);
      n == null || n.search([t, o]).forEach((u) => r.next(u)), r.complete();
    }, s.signal);
  }
  freeResources() {
  }
}
E.capabilities = ["getFeatures", "getRefNames"];
export {
  E as default,
  w as featureData
};
