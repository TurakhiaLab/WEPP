import { aV as n, aZ as l } from "./index-8ee3baa4.js";
import "react";
import "react-dom";
class u extends n.BaseAdapter {
  async getRefNameAliases() {
    const e = this.getConf("location");
    if (e.uri === "" || e.uri === "/path/to/my/aliases.txt")
      return [];
    const s = await l.openLocation(e, this.pluginManager).readFile("utf8"), r = this.getConf("refNameColumn");
    return s.trim().split(/\n|\r\n|\r/).filter((t) => !!t && !t.startsWith("#")).map((t) => {
      const a = t.split("	"), [i] = a.splice(r, 1);
      return { refName: i, aliases: a.filter((o) => !!o.trim()) };
    });
  }
  async freeResources() {
  }
}
export {
  u as default
};
