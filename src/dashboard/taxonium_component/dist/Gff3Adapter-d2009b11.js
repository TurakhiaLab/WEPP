import { aV as h, aZ as p, b9 as m, ar as g, aW as D } from "./index-8ee3baa4.js";
import { I as w } from "./main.esm-53240cd4.js";
import { g as y } from "./index-205eaa5a.js";
import "react";
import "react-dom";
import "./index-aa7a8a76.js";
function b(c) {
  return c[0] === 31 && c[1] === 139 && c[2] === 8;
}
class q extends h.BaseFeatureDataAdapter {
  async loadDataP() {
    const t = this.pluginManager, e = await p.openLocation(this.getConf("gffLocation"), t).readFile(), o = b(e) ? await m(e) : e;
    if (o.length > 536870888)
      throw new Error("Data exceeds maximum string length (512MB)");
    const s = new TextDecoder("utf8", { fatal: !0 }).decode(o), a = s.split(/\n|\r\n|\r/), r = [];
    for (let i = 0; i < a.length && a[i].startsWith("#"); i++)
      r.push(a[i]);
    const n = r.join(`
`), d = y.parseStringSync(s, {
      parseFeatures: !0,
      parseComments: !1,
      parseDirectives: !1,
      parseSequences: !1,
      disableDerivesFromReferences: !0
    }), f = {};
    for (const i of d.flat().map((l, u) => new g({
      data: this.featureData(l),
      id: `${this.id}-offset-${u}`
    }))) {
      const l = i.get("refName");
      f[l] || (f[l] = new w()), f[l].insert([i.get("start"), i.get("end")], i);
    }
    return { header: n, intervalTree: f };
  }
  async loadData() {
    return this.gffFeatures || (this.gffFeatures = this.loadDataP().catch((t) => {
      throw this.gffFeatures = void 0, t;
    })), this.gffFeatures;
  }
  async getRefNames(t = {}) {
    const { intervalTree: e } = await this.loadData();
    return Object.keys(e);
  }
  async getHeader() {
    const { header: t } = await this.loadData();
    return t;
  }
  getFeatures(t, e = {}) {
    return D(async (o) => {
      var s;
      try {
        const { start: a, end: r, refName: n } = t, { intervalTree: d } = await this.loadData();
        (s = d[n]) === null || s === void 0 || s.search([a, r]).forEach((f) => o.next(f)), o.complete();
      } catch (a) {
        o.error(a);
      }
    }, e.signal);
  }
  featureData(t) {
    const e = { ...t };
    e.start -= 1, t.strand === "+" ? e.strand = 1 : t.strand === "-" ? e.strand = -1 : t.strand === "." ? e.strand = 0 : e.strand = void 0, e.phase = Number(t.phase), e.refName = t.seq_id, t.score === null && delete e.score, t.phase === null && delete e.score;
    const o = /* @__PURE__ */ new Set([
      "start",
      "end",
      "seq_id",
      "score",
      "type",
      "source",
      "phase",
      "strand"
    ]), s = t.attributes || {};
    for (const a of Object.keys(s)) {
      let r = a.toLowerCase();
      if (o.has(r) && (r += "2"), s[a] !== null) {
        let n = s[a];
        Array.isArray(n) && n.length === 1 && ([n] = n), e[r] = n;
      }
    }
    return e.refName = e.seq_id, t.child_features && t.child_features.length > 0 && (e.subfeatures = t.child_features.flatMap((a) => a.map((r) => this.featureData(r)))), delete e.child_features, delete e.data, delete e.attributes, delete e.seq_id, e;
  }
  freeResources() {
  }
}
export {
  q as default
};
