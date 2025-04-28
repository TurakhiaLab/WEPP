import { aV as y, bh as p, aZ as x, b9 as N, aW as V } from "./index-8ee3baa4.js";
import { I as b } from "./main.esm-53240cd4.js";
import { V as T } from "./index-7d186318.js";
import "react";
import "react-dom";
const C = (s) => {
  const t = [], a = [];
  return s.split(/\n|\r\n|\r/).map((e) => e.trim()).filter((e) => !!e).forEach((e) => {
    e.startsWith("#") ? t.push(e) : e && a.push(e);
  }), { header: t.join(`
`), lines: a };
};
function D(s) {
  return s[0] === 31 && s[1] === 139 && s[2] === 8;
}
class E extends y.BaseFeatureDataAdapter {
  async getHeader() {
    const { header: t } = await this.setup();
    return t;
  }
  async getMetadata() {
    const { header: t } = await this.setup();
    return new p({ header: t }).getMetadata();
  }
  // converts lines into an interval tree
  async setupP() {
    const t = this.pluginManager, a = await x.openLocation(this.getConf("vcfLocation"), t).readFile(), e = D(a) ? await N(a) : a;
    if (e.length > 536870888)
      throw new Error("Data exceeds maximum string length (512MB)");
    const c = new TextDecoder().decode(e), { header: d, lines: h } = C(c), n = {};
    for (const i of h.map((r, u) => {
      var o;
      const [m, l, , g, , , , v] = r.split("	"), f = +l - 1, w = f + g.length, F = +(((o = v.match(/END=(\d+)/)) === null || o === void 0 ? void 0 : o[1].trim()) || w);
      return { line: r, refName: m, start: f, end: F, id: u };
    })) {
      const r = i.refName;
      n[r] || (n[r] = new b()), n[r].insert([i.start, i.end], i);
    }
    return { header: d, intervalTree: n };
  }
  async setup() {
    return this.vcfFeatures || (this.vcfFeatures = this.setupP().catch((t) => {
      throw this.vcfFeatures = void 0, t;
    })), this.vcfFeatures;
  }
  async getRefNames(t = {}) {
    const { intervalTree: a } = await this.setup();
    return Object.keys(a);
  }
  getFeatures(t, a = {}) {
    return V(async (e) => {
      var c;
      try {
        const { start: d, end: h, refName: n } = t, { header: i, intervalTree: r } = await this.setup(), u = new p({ header: i });
        (c = r[n]) === null || c === void 0 || c.search([d, h]).forEach((o) => e.next(new T({
          variant: u.parseLine(o.line),
          parser: u,
          id: `${this.id}-${o.id}`
        }))), e.complete();
      } catch (d) {
        e.error(d);
      }
    }, a.signal);
  }
  freeResources() {
  }
}
E.capabilities = ["getFeatures", "getRefNames"];
export {
  E as default
};
