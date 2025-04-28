import { Q as c, $ as f } from "./index-8ee3baa4.js";
function g(t, s) {
  for (var o = 0; o < s.length; o++) {
    const e = s[o];
    if (typeof e != "string" && !Array.isArray(e)) {
      for (const r in e)
        if (r !== "default" && !(r in t)) {
          const a = Object.getOwnPropertyDescriptor(e, r);
          a && Object.defineProperty(t, r, a.get ? a : {
            enumerable: !0,
            get: () => e[r]
          });
        }
    }
  }
  return Object.freeze(Object.defineProperty(t, Symbol.toStringTag, { value: "Module" }));
}
var n = f();
const i = /* @__PURE__ */ c(n), u = /* @__PURE__ */ g({
  __proto__: null,
  default: i
}, [n]);
export {
  u as c
};
