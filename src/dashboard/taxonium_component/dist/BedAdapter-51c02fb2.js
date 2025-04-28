import { B as F, f as C } from "./util-1d4160d7.js";
import { aV as y, aZ as x, b9 as I, aW as B } from "./index-8ee3baa4.js";
import { I as E } from "./main.esm-53240cd4.js";
import "react";
import "react-dom";
function R(u) {
  return u[0] === 31 && u[1] === 139 && u[2] === 8;
}
class S extends y.BaseFeatureDataAdapter {
  constructor() {
    super(...arguments), this.intervalTrees = {};
  }
  async loadDataP(e = {}) {
    const t = this.pluginManager, r = this.getConf("bedLocation"), a = await x.openLocation(r, t).readFile(e), n = R(a) ? await I(a) : a;
    if (n.length > 536870888)
      throw new Error("Data exceeds maximum string length (512MB)");
    const s = new TextDecoder("utf8", { fatal: !0 }).decode(n).split(/\n|\r\n|\r/).filter((d) => !!d), i = [];
    let o = 0;
    for (; o < s.length && s[o].startsWith("#"); o++)
      i.push(s[o]);
    const f = i.join(`
`), l = {};
    for (; o < s.length; o++) {
      const d = s[o], D = d.indexOf("	"), g = d.slice(0, D);
      l[g] || (l[g] = []), l[g].push(d);
    }
    const c = this.getConf("autoSql"), m = new F({ autoSql: c }), p = this.getConf("columnNames"), v = this.getConf("scoreColumn"), w = this.getConf("colRef"), T = this.getConf("colStart"), b = this.getConf("colEnd");
    return {
      header: f,
      features: l,
      parser: m,
      columnNames: p,
      scoreColumn: v,
      colRef: w,
      colStart: T,
      colEnd: b
    };
  }
  async loadData(e = {}) {
    return this.bedFeatures || (this.bedFeatures = this.loadDataP(e).catch((t) => {
      throw this.bedFeatures = void 0, t;
    })), this.bedFeatures;
  }
  async getRefNames(e = {}) {
    const { features: t } = await this.loadData(e);
    return Object.keys(t);
  }
  async getHeader(e = {}) {
    const { header: t } = await this.loadData(e);
    return t;
  }
  async getNames() {
    const { header: e, columnNames: t } = await this.loadData();
    if (t.length)
      return t;
    const r = e.split(/\n|\r\n|\r/).filter((n) => !!n), a = r[r.length - 1];
    return a != null && a.includes("	") ? a.slice(1).split("	").map((n) => n.trim()) : void 0;
  }
  async loadFeatureIntervalTreeHelper(e) {
    const { colRef: t, colStart: r, colEnd: a, features: n, parser: h, scoreColumn: s } = await this.loadData(), i = n[e];
    if (!i)
      return;
    const o = await this.getNames(), f = new E(), l = i.map((c, m) => {
      const p = `${this.id}-${e}-${m}`;
      return C(c, t, r, a, s, h, p, o);
    });
    for (const c of l)
      f.insert([c.get("start"), c.get("end")], c);
    return f;
  }
  async loadFeatureIntervalTree(e) {
    return this.intervalTrees[e] || (this.intervalTrees[e] = this.loadFeatureIntervalTreeHelper(e).catch((t) => {
      throw this.intervalTrees[e] = void 0, t;
    })), this.intervalTrees[e];
  }
  getFeatures(e, t = {}) {
    return B(async (r) => {
      const { start: a, end: n, refName: h } = e, s = await this.loadFeatureIntervalTree(h);
      s == null || s.search([a, n]).forEach((i) => r.next(i)), r.complete();
    }, t.signal);
  }
  freeResources() {
  }
}
S.capabilities = ["getFeatures", "getRefNames"];
export {
  S as default
};
