import { bf as c, aZ as s } from "./index-8ee3baa4.js";
import { I as p, a as g } from "./IndexedFastaAdapter-ebec6324.js";
import "react";
import "react-dom";
class d extends p {
  constructor({ fasta: a, path: t, fai: n, faiPath: o, gzi: e, gziPath: i, chunkSizeLimit: r }) {
    super({ fasta: a, path: t, fai: n, faiPath: o, chunkSizeLimit: r }), a && e ? this.fasta = new c({
      filehandle: a,
      gziFilehandle: e
    }) : t && i && (this.fasta = new c({ path: t, gziPath: i }));
  }
}
class m extends g {
  async setupPre() {
    const a = this.getConf("fastaLocation"), t = this.getConf("faiLocation"), n = this.getConf("gziLocation"), o = {
      fasta: s.openLocation(a, this.pluginManager),
      fai: s.openLocation(t, this.pluginManager),
      gzi: s.openLocation(n, this.pluginManager)
    };
    return { fasta: new d(o) };
  }
}
export {
  m as default
};
