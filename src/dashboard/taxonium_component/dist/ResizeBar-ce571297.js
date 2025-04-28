import { a1 as c, b2 as b, b0 as p } from "./index-8ee3baa4.js";
import g from "react";
var f = {}, m = c && c.__createBinding || (Object.create ? function(e, t, r, n) {
  n === void 0 && (n = r);
  var a = Object.getOwnPropertyDescriptor(t, r);
  (!a || ("get" in a ? !t.__esModule : a.writable || a.configurable)) && (a = { enumerable: !0, get: function() {
    return t[r];
  } }), Object.defineProperty(e, n, a);
} : function(e, t, r, n) {
  n === void 0 && (n = r), e[n] = t[r];
}), h = c && c.__setModuleDefault || (Object.create ? function(e, t) {
  Object.defineProperty(e, "default", { enumerable: !0, value: t });
} : function(e, t) {
  e.default = t;
}), y = c && c.__importStar || function(e) {
  if (e && e.__esModule)
    return e;
  var t = {};
  if (e != null)
    for (var r in e)
      r !== "default" && Object.prototype.hasOwnProperty.call(e, r) && m(t, e, r);
  return h(t, e), t;
}, B = c && c.__importDefault || function(e) {
  return e && e.__esModule ? e : { default: e };
};
Object.defineProperty(f, "__esModule", { value: !0 });
var z = f.useResizeBar = void 0;
const l = y(g), R = p, O = B(b), _ = (0, R.makeStyles)()((e) => ({
  resizeBar: {
    background: e.palette.action.disabledBackground,
    height: 12,
    position: "relative",
    overflow: "hidden"
  },
  tick: {
    position: "absolute",
    height: "100%",
    pointerEvents: "none",
    background: e.palette.action.disabled,
    width: 1
  },
  hiddenTick: {
    position: "absolute",
    height: "100%",
    width: 5
  }
}));
function j() {
  const e = (0, l.useRef)(null), [t, r] = (0, l.useState)(0);
  return (0, l.useEffect)(() => {
    const n = setInterval(() => {
      var a;
      const s = (a = e.current) === null || a === void 0 ? void 0 : a.querySelector(".MuiDataGrid-virtualScroller");
      s && r(s.scrollLeft);
    }, 100);
    return () => {
      clearInterval(n);
    };
  }, []), { ref: e, scrollLeft: t };
}
z = f.useResizeBar = j;
function D({ left: e, scrollLeft: t, idx: r, onDrag: n }) {
  const { classes: a } = _(), s = (0, l.useCallback)((o) => {
    n(o, r);
  }, [r, n]);
  return l.default.createElement(
    l.default.Fragment,
    null,
    l.default.createElement(O.default, { onDrag: s, vertical: !0, className: a.hiddenTick, style: { left: e - t - 2.5 } }),
    l.default.createElement("div", { style: { left: e - t }, className: a.tick })
  );
}
function M({ widths: e, setWidths: t, checkbox: r, scrollLeft: n = 0 }) {
  const { classes: a } = _(), s = [];
  let o = r ? 52 : 0;
  for (let i = 0; i < e.length; i++) {
    const u = e[i];
    s[i] = u + o, o += u;
  }
  const v = (0, l.useCallback)((i, u) => {
    const d = [...e];
    d[u] = Math.max(d[u] + i, 50), t(d);
  }, [e, t]);
  return l.default.createElement("div", { className: a.resizeBar }, s.map((i, u) => l.default.createElement(D, { key: u, left: u === s.length - 1 ? i - 3 : i, onDrag: v, idx: u, scrollLeft: n })));
}
var w = f.default = M;
export {
  w as _,
  z as u
};
