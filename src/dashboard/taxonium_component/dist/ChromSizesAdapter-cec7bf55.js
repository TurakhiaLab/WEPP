import { aV as n, aZ as o } from "./index-8ee3baa4.js";
import "react";
import "react-dom";
class m extends n.BaseAdapter {
  async setupPre() {
    const t = this.pluginManager, r = await o.openLocation(this.getConf("chromSizesLocation"), t).readFile("utf8");
    return Object.fromEntries(r.split(/\n|\r\n|\r/).map((e) => e.trim()).filter((e) => !!e).map((e) => {
      const [a, i] = e.split("	");
      return [a, +i];
    }));
  }
  async setup() {
    return this.setupP || (this.setupP = this.setupPre().catch((t) => {
      throw this.setupP = void 0, t;
    })), this.setupP;
  }
  async getRegions() {
    const t = await this.setup();
    return Object.keys(t).map((s) => ({
      refName: s,
      start: 0,
      end: t[s]
    }));
  }
  getHeader() {
    return {};
  }
  freeResources() {
  }
}
export {
  m as default
};
