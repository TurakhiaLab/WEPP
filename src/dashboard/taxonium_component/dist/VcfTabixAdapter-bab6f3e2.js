import { aV as u, aZ as o, bh as g, aW as h } from "./index-8ee3baa4.js";
import { T as p } from "./tabixIndexedFile-24223997.js";
import "./long-a46ae85c.js";
import { V as l } from "./index-7d186318.js";
import "react";
import "react-dom";
import "./index-05a6222b.js";
class b extends u.BaseFeatureDataAdapter {
  async configurePre() {
    const e = this.pluginManager, t = this.getConf("vcfGzLocation"), a = this.getConf(["index", "location"]), s = this.getConf(["index", "indexType"]), c = o.openLocation(t, e), n = s === "CSI", i = new p({
      filehandle: c,
      csiFilehandle: n ? o.openLocation(a, e) : void 0,
      tbiFilehandle: n ? void 0 : o.openLocation(a, e),
      chunkCacheSize: 50 * 2 ** 20,
      chunkSizeLimit: 1e9
    }), r = await i.getHeader();
    return {
      vcf: i,
      parser: new g({ header: r })
    };
  }
  async configure() {
    return this.configured || (this.configured = this.configurePre().catch((e) => {
      throw this.configured = void 0, e;
    })), this.configured;
  }
  async getRefNames(e = {}) {
    const { vcf: t } = await this.configure();
    return t.getReferenceSequenceNames(e);
  }
  async getHeader() {
    const { vcf: e } = await this.configure();
    return e.getHeader();
  }
  async getMetadata() {
    const { parser: e } = await this.configure();
    return e.getMetadata();
  }
  getFeatures(e, t = {}) {
    return h(async (a) => {
      const { refName: s, start: c, end: n } = e, { vcf: i, parser: r } = await this.configure();
      await i.getLines(s, c, n, {
        lineCallback: (f, d) => {
          a.next(new l({
            variant: r.parseLine(f),
            parser: r,
            id: `${this.id}-vcf-${d}`
          }));
        },
        ...t
      }), a.complete();
    }, t.signal);
  }
  freeResources() {
  }
}
export {
  b as default
};
