import { aV as g, l, aZ as h, aW as m, be as b, ar as _ } from "./index-8ee3baa4.js";
import { T as F } from "./tabixIndexedFile-24223997.js";
import "./long-a46ae85c.js";
import { g as H } from "./index-205eaa5a.js";
import "react";
import "react-dom";
import "./index-05a6222b.js";
import "./index-aa7a8a76.js";
class N extends g.BaseFeatureDataAdapter {
  constructor(e, t, r) {
    super(e, t, r);
    const a = l.readConfObject(e, "gffGzLocation"), f = l.readConfObject(e, ["index", "indexType"]), i = l.readConfObject(e, ["index", "location"]), n = l.readConfObject(e, "dontRedispatch");
    this.dontRedispatch = n || ["chromosome", "contig", "region"], this.gff = new F({
      filehandle: h.openLocation(a, this.pluginManager),
      csiFilehandle: f === "CSI" ? h.openLocation(i, this.pluginManager) : void 0,
      tbiFilehandle: f !== "CSI" ? h.openLocation(i, this.pluginManager) : void 0,
      chunkCacheSize: 50 * 2 ** 20,
      renameRefSeqs: (c) => c
    });
  }
  async getRefNames(e = {}) {
    return this.gff.getReferenceSequenceNames(e);
  }
  async getHeader() {
    return this.gff.getHeader();
  }
  getFeatures(e, t = {}) {
    return m(async (r) => {
      const a = await this.gff.getMetadata();
      await this.getFeaturesHelper(e, t, a, r, !0);
    }, t.signal);
  }
  async getFeaturesHelper(e, t = {}, r, a, f, i = e) {
    try {
      const n = [];
      if (await this.gff.getLines(e.refName, e.start, e.end, (s, d) => {
        n.push(this.parseLine(r.columnNumbers, s, d));
      }), f && n.length) {
        let s = 1 / 0, d = -1 / 0;
        if (n.forEach((o) => {
          const p = o.fields[2];
          if (!this.dontRedispatch.includes(p)) {
            const u = o.start - 1;
            u < s && (s = u), o.end > d && (d = o.end);
          }
        }), d > e.end || s < e.start) {
          await this.getFeaturesHelper({ ...e, start: s, end: d }, t, r, a, !1, e);
          return;
        }
      }
      const c = n.map((s) => (s.fields[8] && s.fields[8] !== "." ? s.fields[8].includes("_lineHash") || (s.fields[8] += `;_lineHash=${s.lineHash}`) : s.fields[8] = `_lineHash=${s.lineHash}`, s.fields.join("	"))).join(`
`);
      H.parseStringSync(c, {
        parseFeatures: !0,
        parseComments: !1,
        parseDirectives: !1,
        parseSequences: !1,
        disableDerivesFromReferences: !0
      }).forEach((s) => this.formatFeatures(s).forEach((d) => {
        b(d.get("start"), d.get("end"), i.start, i.end) && a.next(d);
      })), a.complete();
    } catch (n) {
      a.error(n);
    }
  }
  parseLine(e, t, r) {
    const a = t.split("	");
    return {
      start: +a[e.start - 1],
      end: +a[e.end - 1],
      lineHash: r,
      fields: a
    };
  }
  formatFeatures(e) {
    return e.map((t) => new _({
      data: this.featureData(t),
      // eslint-disable-next-line @typescript-eslint/no-non-null-assertion
      id: `${this.id}-offset-${t.attributes._lineHash[0]}`
    }));
  }
  featureData(e) {
    const t = { ...e };
    t.start -= 1, e.strand === "+" ? t.strand = 1 : e.strand === "-" ? t.strand = -1 : e.strand === "." ? t.strand = 0 : t.strand = void 0, t.phase = Number(e.phase), t.refName = e.seq_id, e.score === null && delete t.score, e.phase === null && delete t.score;
    const r = /* @__PURE__ */ new Set([
      "start",
      "end",
      "seq_id",
      "score",
      "type",
      "source",
      "phase",
      "strand"
    ]), a = e.attributes || {};
    for (const f of Object.keys(a)) {
      let i = f.toLowerCase();
      if (r.has(i) && (i += "2"), a[f] !== null) {
        let n = a[f];
        Array.isArray(n) && n.length === 1 && ([n] = n), t[i] = n;
      }
    }
    return t.refName = t.seq_id, e.child_features && e.child_features.length > 0 && (t.subfeatures = e.child_features.flatMap((f) => f.map((i) => this.featureData(i)))), delete t.child_features, delete t.data, delete t._linehash, delete t.attributes, delete t.seq_id, t;
  }
  freeResources() {
  }
}
export {
  N as default
};
