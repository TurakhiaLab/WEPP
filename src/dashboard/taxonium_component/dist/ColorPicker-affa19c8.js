import { a1 as C, a_ as pe, a$ as me, Q as be, b0 as _e, b1 as Ce } from "./index-8ee3baa4.js";
import re from "react";
import { c as Ee } from "./index-090e001f.js";
function Pe(e, t) {
  for (var r = 0; r < t.length; r++) {
    const o = t[r];
    if (typeof o != "string" && !Array.isArray(o)) {
      for (const a in o)
        if (a !== "default" && !(a in e)) {
          const l = Object.getOwnPropertyDescriptor(o, a);
          l && Object.defineProperty(e, a, l.get ? l : {
            enumerable: !0,
            get: () => o[a]
          });
        }
    }
  }
  return Object.freeze(Object.defineProperty(e, Symbol.toStringTag, { value: "Module" }));
}
var G = {}, i = {}, ye = C && C.__createBinding || (Object.create ? function(e, t, r, o) {
  o === void 0 && (o = r);
  var a = Object.getOwnPropertyDescriptor(t, r);
  (!a || ("get" in a ? !t.__esModule : a.writable || a.configurable)) && (a = { enumerable: !0, get: function() {
    return t[r];
  } }), Object.defineProperty(e, o, a);
} : function(e, t, r, o) {
  o === void 0 && (o = r), e[o] = t[r];
}), He = C && C.__setModuleDefault || (Object.create ? function(e, t) {
  Object.defineProperty(e, "default", { enumerable: !0, value: t });
} : function(e, t) {
  e.default = t;
}), ke = C && C.__importStar || function(e) {
  if (e && e.__esModule)
    return e;
  var t = {};
  if (e != null)
    for (var r in e)
      r !== "default" && Object.prototype.hasOwnProperty.call(e, r) && ye(t, e, r);
  return He(t, e), t;
};
Object.defineProperty(i, "__esModule", { value: !0 });
i.setNonce = i.RgbaStringColorPicker = i.RgbaColorPicker = i.RgbStringColorPicker = i.RgbColorPicker = i.HsvaStringColorPicker = i.HsvaColorPicker = i.HsvStringColorPicker = i.HsvColorPicker = i.HslaStringColorPicker = i.HslaColorPicker = i.HslStringColorPicker = i.HslColorPicker = i.HexColorPicker = i.HexColorInput = void 0;
const n = ke(re);
function p() {
  return (p = Object.assign || function(e) {
    for (var t = 1; t < arguments.length; t++) {
      var r = arguments[t];
      for (var o in r)
        Object.prototype.hasOwnProperty.call(r, o) && (e[o] = r[o]);
    }
    return e;
  }).apply(this, arguments);
}
function R(e, t) {
  if (e == null)
    return {};
  var r, o, a = {}, l = Object.keys(e);
  for (o = 0; o < l.length; o++)
    t.indexOf(r = l[o]) >= 0 || (a[r] = e[r]);
  return a;
}
function j(e) {
  const t = (0, n.useRef)(e), r = (0, n.useRef)((o) => {
    t.current && t.current(o);
  });
  return t.current = e, r.current;
}
const w = (e, t = 0, r = 1) => e > r ? r : e < t ? t : e, S = (e) => "touches" in e, Y = (e) => e && e.ownerDocument.defaultView || self, J = (e, t, r) => {
  const o = e.getBoundingClientRect(), a = S(t) ? ((l, s) => {
    for (let u = 0; u < l.length; u++)
      if (l[u].identifier === s)
        return l[u];
    return l[0];
  })(t.touches, r) : t;
  return {
    left: w((a.pageX - (o.left + Y(e).pageXOffset)) / o.width),
    top: w((a.pageY - (o.top + Y(e).pageYOffset)) / o.height)
  };
}, Z = (e) => {
  !S(e) && e.preventDefault();
}, Q = n.default.memo((e) => {
  let { onMove: t, onKey: r } = e, o = R(e, ["onMove", "onKey"]);
  const a = (0, n.useRef)(null), l = j(t), s = j(r), u = (0, n.useRef)(null), f = (0, n.useRef)(!1), [m, E, k] = (0, n.useMemo)(() => {
    const M = (g) => {
      Z(g), (S(g) ? g.touches.length > 0 : g.buttons > 0) && a.current ? l(J(a.current, g, u.current)) : y(!1);
    }, O = () => y(!1);
    function y(g) {
      const c = f.current, d = Y(a.current), h = g ? d.addEventListener : d.removeEventListener;
      h(c ? "touchmove" : "mousemove", M), h(c ? "touchend" : "mouseup", O);
    }
    return [
      ({ nativeEvent: g }) => {
        const c = a.current;
        if (c && (Z(g), !((d, h) => h && !S(d))(g, f.current) && c)) {
          if (S(g)) {
            f.current = !0;
            const d = g.changedTouches || [];
            d.length && (u.current = d[0].identifier);
          }
          c.focus(), l(J(c, g, u.current)), y(!0);
        }
      },
      (g) => {
        const c = g.which || g.keyCode;
        c < 37 || c > 40 || (g.preventDefault(), s({
          left: c === 39 ? 0.05 : c === 37 ? -0.05 : 0,
          top: c === 40 ? 0.05 : c === 38 ? -0.05 : 0
        }));
      },
      y
    ];
  }, [s, l]);
  return (0, n.useEffect)(() => k, [k]), n.default.createElement("div", p({}, o, {
    onTouchStart: m,
    onMouseDown: m,
    className: "react-colorful__interactive",
    ref: a,
    onKeyDown: E,
    tabIndex: 0,
    role: "slider"
  }));
}), q = (e) => e.filter(Boolean).join(" "), U = ({ className: e, color: t, left: r, top: o = 0.5 }) => {
  const a = q(["react-colorful__pointer", e]);
  return n.default.createElement("div", { className: a, style: { top: 100 * o + "%", left: 100 * r + "%" } }, n.default.createElement("div", {
    className: "react-colorful__pointer-fill",
    style: { backgroundColor: t }
  }));
}, v = (e, t = 0, r = Math.pow(10, t)) => Math.round(r * e) / r, Me = { grad: 0.9, turn: 360, rad: 360 / (2 * Math.PI) }, V = (e) => (e[0] === "#" && (e = e.substr(1)), e.length < 6 ? {
  r: parseInt(e[0] + e[0], 16),
  g: parseInt(e[1] + e[1], 16),
  b: parseInt(e[2] + e[2], 16),
  a: 1
} : {
  r: parseInt(e.substr(0, 2), 16),
  g: parseInt(e.substr(2, 2), 16),
  b: parseInt(e.substr(4, 2), 16),
  a: 1
}), oe = (e, t = "deg") => Number(e) * (Me[t] || 1), ae = (e) => {
  const t = /hsla?\(?\s*(-?\d*\.?\d+)(deg|rad|grad|turn)?[,\s]+(-?\d*\.?\d+)%?[,\s]+(-?\d*\.?\d+)%?,?\s*[/\s]*(-?\d*\.?\d+)?(%)?\s*\)?/i.exec(e);
  return t ? W({
    h: oe(t[1], t[2]),
    s: Number(t[3]),
    l: Number(t[4]),
    a: t[5] === void 0 ? 1 : Number(t[5]) / (t[6] ? 100 : 1)
  }) : { h: 0, s: 0, v: 0, a: 1 };
}, xe = ae, W = ({ h: e, s: t, l: r, a: o }) => ({
  h: e,
  s: (t *= (r < 50 ? r : 100 - r) / 100) > 0 ? 2 * t / (r + t) * 100 : 0,
  v: r + t,
  a: o
}), A = ({ h: e, s: t, v: r, a: o }) => {
  const a = (200 - t) * r / 100;
  return {
    h: v(e),
    s: v(a > 0 && a < 200 ? t * r / 100 / (a <= 100 ? a : 200 - a) * 100 : 0),
    l: v(a / 2),
    a: v(o, 2)
  };
}, T = (e) => {
  const { h: t, s: r, l: o } = A(e);
  return `hsl(${t}, ${r}%, ${o}%)`;
}, z = (e) => {
  const { h: t, s: r, l: o, a } = A(e);
  return `hsla(${t}, ${r}%, ${o}%, ${a})`;
}, D = ({ h: e, s: t, v: r, a: o }) => {
  e = e / 360 * 6, t /= 100, r /= 100;
  const a = Math.floor(e), l = r * (1 - t), s = r * (1 - (e - a) * t), u = r * (1 - (1 - e + a) * t), f = a % 6;
  return {
    r: v(255 * [r, s, l, l, u, r][f]),
    g: v(255 * [u, r, r, s, l, l][f]),
    b: v(255 * [l, l, u, r, r, s][f]),
    a: v(o, 2)
  };
}, le = (e) => {
  const t = /hsva?\(?\s*(-?\d*\.?\d+)(deg|rad|grad|turn)?[,\s]+(-?\d*\.?\d+)%?[,\s]+(-?\d*\.?\d+)%?,?\s*[/\s]*(-?\d*\.?\d+)?(%)?\s*\)?/i.exec(e);
  return t ? I({
    h: oe(t[1], t[2]),
    s: Number(t[3]),
    v: Number(t[4]),
    a: t[5] === void 0 ? 1 : Number(t[5]) / (t[6] ? 100 : 1)
  }) : { h: 0, s: 0, v: 0, a: 1 };
}, we = le, ne = (e) => {
  const t = /rgba?\(?\s*(-?\d*\.?\d+)(%)?[,\s]+(-?\d*\.?\d+)(%)?[,\s]+(-?\d*\.?\d+)(%)?,?\s*[/\s]*(-?\d*\.?\d+)?(%)?\s*\)?/i.exec(e);
  return t ? L({
    r: Number(t[1]) / (t[2] ? 100 / 255 : 1),
    g: Number(t[3]) / (t[4] ? 100 / 255 : 1),
    b: Number(t[5]) / (t[6] ? 100 / 255 : 1),
    a: t[7] === void 0 ? 1 : Number(t[7]) / (t[8] ? 100 : 1)
  }) : { h: 0, s: 0, v: 0, a: 1 };
}, Ne = ne, X = (e) => {
  const t = e.toString(16);
  return t.length < 2 ? "0" + t : t;
}, L = ({ r: e, g: t, b: r, a: o }) => {
  const a = Math.max(e, t, r), l = a - Math.min(e, t, r), s = l ? a === e ? (t - r) / l : a === t ? 2 + (r - e) / l : 4 + (e - t) / l : 0;
  return {
    h: v(60 * (s < 0 ? s + 6 : s)),
    s: v(a ? l / a * 100 : 0),
    v: v(a / 255 * 100),
    a: o
  };
}, I = (e) => ({ h: v(e.h), s: v(e.s), v: v(e.v), a: v(e.a, 2) }), ce = n.default.memo(({ className: e, hue: t, onChange: r }) => {
  const o = q(["react-colorful__hue", e]);
  return n.default.createElement("div", { className: o }, n.default.createElement(Q, {
    onMove: (a) => {
      r({ h: 360 * a.left });
    },
    onKey: (a) => {
      r({ h: w(t + 360 * a.left, 0, 360) });
    },
    "aria-label": "Hue",
    "aria-valuetext": v(t)
  }, n.default.createElement(U, {
    className: "react-colorful__hue-pointer",
    left: t / 360,
    color: T({ h: t, s: 100, v: 100, a: 1 })
  })));
}), se = n.default.memo(({ hsva: e, onChange: t }) => {
  const r = { backgroundColor: T({ h: e.h, s: 100, v: 100, a: 1 }) };
  return n.default.createElement("div", { className: "react-colorful__saturation", style: r }, n.default.createElement(Q, {
    onMove: (o) => {
      t({ s: 100 * o.left, v: 100 - 100 * o.top });
    },
    onKey: (o) => {
      t({
        s: w(e.s + 100 * o.left, 0, 100),
        v: w(e.v - 100 * o.top, 0, 100)
      });
    },
    "aria-label": "Color",
    "aria-valuetext": `Saturation ${v(e.s)}%, Brightness ${v(e.v)}%`
  }, n.default.createElement(U, {
    className: "react-colorful__saturation-pointer",
    top: 1 - e.v / 100,
    left: e.s / 100,
    color: T(e)
  })));
}), H = (e, t) => {
  if (e === t)
    return !0;
  for (const r in e)
    if (e[r] !== t[r])
      return !1;
  return !0;
}, N = (e, t) => e.replace(/\s/g, "") === t.replace(/\s/g, "");
function ue(e, t, r) {
  const o = j(r), [a, l] = (0, n.useState)(() => e.toHsva(t)), s = (0, n.useRef)({ color: t, hsva: a });
  (0, n.useEffect)(() => {
    if (!e.equal(t, s.current.color)) {
      const f = e.toHsva(t);
      s.current = { hsva: f, color: t }, l(f);
    }
  }, [t, e]), (0, n.useEffect)(() => {
    let f;
    H(a, s.current.hsva) || e.equal(f = e.fromHsva(a), s.current.color) || (s.current = { hsva: a, color: f }, o(f));
  }, [a, e, o]);
  const u = (0, n.useCallback)((f) => {
    l((m) => Object.assign({}, m, f));
  }, []);
  return [a, u];
}
const $e = typeof window < "u" ? n.useLayoutEffect : n.useEffect;
let ie;
const Oe = () => ie || (typeof __webpack_nonce__ < "u" ? __webpack_nonce__ : void 0), Se = (e) => {
  ie = e;
}, ee = /* @__PURE__ */ new Map(), fe = (e) => {
  $e(() => {
    const t = e.current ? e.current.ownerDocument : document;
    if (t !== void 0 && !ee.has(t)) {
      const r = t.createElement("style");
      r.innerHTML = `.react-colorful{position:relative;display:flex;flex-direction:column;width:200px;height:200px;-webkit-user-select:none;-moz-user-select:none;-ms-user-select:none;user-select:none;cursor:default}.react-colorful__saturation{position:relative;flex-grow:1;border-color:transparent;border-bottom:12px solid #000;border-radius:8px 8px 0 0;background-image:linear-gradient(0deg,#000,transparent),linear-gradient(90deg,#fff,hsla(0,0%,100%,0))}.react-colorful__alpha-gradient,.react-colorful__pointer-fill{content:"";position:absolute;left:0;top:0;right:0;bottom:0;pointer-events:none;border-radius:inherit}.react-colorful__alpha-gradient,.react-colorful__saturation{box-shadow:inset 0 0 0 1px rgba(0,0,0,.05)}.react-colorful__alpha,.react-colorful__hue{position:relative;height:24px}.react-colorful__hue{background:linear-gradient(90deg,red 0,#ff0 17%,#0f0 33%,#0ff 50%,#00f 67%,#f0f 83%,red)}.react-colorful__last-control{border-radius:0 0 8px 8px}.react-colorful__interactive{position:absolute;left:0;top:0;right:0;bottom:0;border-radius:inherit;outline:none;touch-action:none}.react-colorful__pointer{position:absolute;z-index:1;box-sizing:border-box;width:28px;height:28px;transform:translate(-50%,-50%);background-color:#fff;border:2px solid #fff;border-radius:50%;box-shadow:0 2px 4px rgba(0,0,0,.2)}.react-colorful__interactive:focus .react-colorful__pointer{transform:translate(-50%,-50%) scale(1.1)}.react-colorful__alpha,.react-colorful__alpha-pointer{background-color:#fff;background-image:url('data:image/svg+xml;charset=utf-8,<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill-opacity=".05"><path d="M8 0h8v8H8zM0 8h8v8H0z"/></svg>')}.react-colorful__saturation-pointer{z-index:3}.react-colorful__hue-pointer{z-index:2}`, ee.set(t, r);
      const o = Oe();
      o && r.setAttribute("nonce", o), t.head.appendChild(r);
    }
  }, []);
}, x = (e) => {
  let { className: t, colorModel: r, color: o = r.defaultColor, onChange: a } = e, l = R(e, ["className", "colorModel", "color", "onChange"]);
  const s = (0, n.useRef)(null);
  fe(s);
  const [u, f] = ue(r, o, a), m = q(["react-colorful", t]);
  return n.default.createElement("div", p({}, l, { ref: s, className: m }), n.default.createElement(se, { hsva: u, onChange: f }), n.default.createElement(ce, {
    hue: u.h,
    onChange: f,
    className: "react-colorful__last-control"
  }));
}, je = {
  defaultColor: "000",
  toHsva: (e) => L(V(e)),
  fromHsva: (e) => (({ r: t, g: r, b: o }) => "#" + X(t) + X(r) + X(o))(D(e)),
  equal: (e, t) => e.toLowerCase() === t.toLowerCase() || H(V(e), V(t))
}, Re = (e) => n.default.createElement(x, p({}, e, { colorModel: je })), qe = ({ className: e, hsva: t, onChange: r }) => {
  const o = {
    backgroundImage: `linear-gradient(90deg, ${z(Object.assign({}, t, { a: 0 }))}, ${z(Object.assign({}, t, { a: 1 }))})`
  }, a = q(["react-colorful__alpha", e]);
  return n.default.createElement("div", { className: a }, n.default.createElement("div", {
    className: "react-colorful__alpha-gradient",
    style: o
  }), n.default.createElement(Q, {
    onMove: (l) => {
      r({ a: l.left });
    },
    onKey: (l) => {
      r({ a: w(t.a + l.left) });
    },
    "aria-label": "Alpha",
    "aria-valuetext": `${v(100 * t.a)}%`
  }, n.default.createElement(U, {
    className: "react-colorful__alpha-pointer",
    left: t.a,
    color: z(t)
  })));
}, $ = (e) => {
  let { className: t, colorModel: r, color: o = r.defaultColor, onChange: a } = e, l = R(e, ["className", "colorModel", "color", "onChange"]);
  const s = (0, n.useRef)(null);
  fe(s);
  const [u, f] = ue(r, o, a), m = q(["react-colorful", t]);
  return n.default.createElement("div", p({}, l, { ref: s, className: m }), n.default.createElement(se, { hsva: u, onChange: f }), n.default.createElement(ce, { hue: u.h, onChange: f }), n.default.createElement(qe, {
    hsva: u,
    onChange: f,
    className: "react-colorful__last-control"
  }));
}, De = {
  defaultColor: { h: 0, s: 0, l: 0, a: 1 },
  toHsva: W,
  fromHsva: A,
  equal: H
}, Ie = (e) => n.default.createElement($, p({}, e, { colorModel: De })), Be = { defaultColor: "hsla(0, 0%, 0%, 1)", toHsva: ae, fromHsva: z, equal: N }, ze = (e) => n.default.createElement($, p({}, e, { colorModel: Be })), Te = {
  defaultColor: { h: 0, s: 0, l: 0 },
  toHsva: ({ h: e, s: t, l: r }) => W({ h: e, s: t, l: r, a: 1 }),
  fromHsva: (e) => (({ h: t, s: r, l: o }) => ({ h: t, s: r, l: o }))(A(e)),
  equal: H
}, Ae = (e) => n.default.createElement(x, p({}, e, { colorModel: Te })), Le = { defaultColor: "hsl(0, 0%, 0%)", toHsva: xe, fromHsva: T, equal: N }, Ke = (e) => n.default.createElement(x, p({}, e, { colorModel: Le })), Fe = {
  defaultColor: { h: 0, s: 0, v: 0, a: 1 },
  toHsva: (e) => e,
  fromHsva: I,
  equal: H
}, Ve = (e) => n.default.createElement($, p({}, e, { colorModel: Fe })), Xe = {
  defaultColor: "hsva(0, 0%, 0%, 1)",
  toHsva: le,
  fromHsva: (e) => {
    const { h: t, s: r, v: o, a } = I(e);
    return `hsva(${t}, ${r}%, ${o}%, ${a})`;
  },
  equal: N
}, Ye = (e) => n.default.createElement($, p({}, e, { colorModel: Xe })), Ge = {
  defaultColor: { h: 0, s: 0, v: 0 },
  toHsva: ({ h: e, s: t, v: r }) => ({ h: e, s: t, v: r, a: 1 }),
  fromHsva: (e) => {
    const { h: t, s: r, v: o } = I(e);
    return { h: t, s: r, v: o };
  },
  equal: H
}, Qe = (e) => n.default.createElement(x, p({}, e, { colorModel: Ge })), Ue = {
  defaultColor: "hsv(0, 0%, 0%)",
  toHsva: we,
  fromHsva: (e) => {
    const { h: t, s: r, v: o } = I(e);
    return `hsv(${t}, ${r}%, ${o}%)`;
  },
  equal: N
}, We = (e) => n.default.createElement(x, p({}, e, { colorModel: Ue })), Je = {
  defaultColor: { r: 0, g: 0, b: 0, a: 1 },
  toHsva: L,
  fromHsva: D,
  equal: H
}, Ze = (e) => n.default.createElement($, p({}, e, { colorModel: Je })), et = {
  defaultColor: "rgba(0, 0, 0, 1)",
  toHsva: ne,
  fromHsva: (e) => {
    const { r: t, g: r, b: o, a } = D(e);
    return `rgba(${t}, ${r}, ${o}, ${a})`;
  },
  equal: N
}, tt = (e) => n.default.createElement($, p({}, e, { colorModel: et })), rt = {
  defaultColor: { r: 0, g: 0, b: 0 },
  toHsva: ({ r: e, g: t, b: r }) => L({ r: e, g: t, b: r, a: 1 }),
  fromHsva: (e) => (({ r: t, g: r, b: o }) => ({ r: t, g: r, b: o }))(D(e)),
  equal: H
}, ot = (e) => n.default.createElement(x, p({}, e, { colorModel: rt })), at = {
  defaultColor: "rgb(0, 0, 0)",
  toHsva: Ne,
  fromHsva: (e) => {
    const { r: t, g: r, b: o } = D(e);
    return `rgb(${t}, ${r}, ${o})`;
  },
  equal: N
}, lt = (e) => n.default.createElement(x, p({}, e, { colorModel: at })), nt = /^#?([0-9A-F]{3,8})$/i, ct = (e) => {
  const { color: t = "", onChange: r, onBlur: o, escape: a, validate: l, format: s, process: u } = e, f = R(e, [
    "color",
    "onChange",
    "onBlur",
    "escape",
    "validate",
    "format",
    "process"
  ]), [m, E] = (0, n.useState)(() => a(t)), k = j(r), M = j(o), O = (0, n.useCallback)((g) => {
    const c = a(g.target.value);
    E(c), l(c) && k(u ? u(c) : c);
  }, [a, u, l, k]), y = (0, n.useCallback)((g) => {
    l(g.target.value) || E(a(t)), M(g);
  }, [t, a, l, M]);
  return (0, n.useEffect)(() => {
    E(a(t));
  }, [t, a]), n.default.createElement("input", p({}, f, {
    value: s ? s(m) : m,
    spellCheck: "false",
    onChange: O,
    onBlur: y
  }));
}, te = (e) => "#" + e, st = (e) => {
  const { prefixed: t, alpha: r } = e, o = R(e, ["prefixed", "alpha"]), a = (0, n.useCallback)((s) => s.replace(/([^0-9A-F]+)/gi, "").substr(0, r ? 8 : 6), [r]), l = (0, n.useCallback)((s) => ((u, f) => {
    const m = nt.exec(u), E = m ? m[1].length : 0;
    return E === 3 || E === 6 || !!f && E === 4 || !!f && E === 8;
  })(s, r), [r]);
  return n.default.createElement(ct, p({}, o, {
    escape: a,
    format: t ? te : void 0,
    process: te,
    validate: l
  }));
};
i.setNonce = Se;
i.HexColorPicker = Re;
i.HslaColorPicker = Ie;
i.HslaStringColorPicker = ze;
i.HslColorPicker = Ae;
i.HslStringColorPicker = Ke;
i.HsvaColorPicker = Ve;
i.HsvaStringColorPicker = Ye;
i.HsvColorPicker = Qe;
i.HsvStringColorPicker = We;
i.RgbaColorPicker = Ze;
i.RgbaStringColorPicker = tt;
i.RgbColorPicker = ot;
i.RgbStringColorPicker = lt;
i.HexColorInput = st;
(function(e) {
  var t = C && C.__createBinding || (Object.create ? function(c, d, h, b) {
    b === void 0 && (b = h);
    var P = Object.getOwnPropertyDescriptor(d, h);
    (!P || ("get" in P ? !d.__esModule : P.writable || P.configurable)) && (P = { enumerable: !0, get: function() {
      return d[h];
    } }), Object.defineProperty(c, b, P);
  } : function(c, d, h, b) {
    b === void 0 && (b = h), c[b] = d[h];
  }), r = C && C.__setModuleDefault || (Object.create ? function(c, d) {
    Object.defineProperty(c, "default", { enumerable: !0, value: d });
  } : function(c, d) {
    c.default = d;
  }), o = C && C.__importStar || function(c) {
    if (c && c.__esModule)
      return c;
    var d = {};
    if (c != null)
      for (var h in c)
        h !== "default" && Object.prototype.hasOwnProperty.call(c, h) && t(d, c, h);
    return r(d, c), d;
  }, a = C && C.__importDefault || function(c) {
    return c && c.__esModule ? c : { default: c };
  };
  Object.defineProperty(e, "__esModule", { value: !0 }), e.ColorPicker = e.ColorPopover = e.PopoverPicker = void 0;
  const l = o(re), s = a(Ee), u = Ce, f = _e, m = o(pe), E = me(), k = i, M = (0, f.makeStyles)()({
    picker: { position: "relative" },
    swatches: {
      display: "flex",
      padding: 12,
      flexWrap: "wrap"
    },
    swatch: {
      width: 24,
      height: 24,
      margin: 4,
      border: "none",
      padding: 0,
      cursor: "pointer",
      outline: "none"
    }
  }), O = ({ color: c, onChange: d }) => {
    const [h, b] = (0, l.useState)(null), { classes: P } = M();
    return l.default.createElement(
      "div",
      { className: P.picker },
      l.default.createElement("div", { className: P.swatch, style: { backgroundColor: c }, onClick: (K) => b(K.currentTarget) }),
      l.default.createElement(y, { anchorEl: h, onClose: () => b(null), color: c, onChange: d })
    );
  };
  e.PopoverPicker = O;
  function y({ anchorEl: c, onChange: d, onClose: h, color: b }) {
    return l.default.createElement(
      u.Popover,
      { open: !!c, anchorEl: c, onClose: h },
      l.default.createElement(g, { color: b, onChange: d })
    );
  }
  e.ColorPopover = y;
  function g({ onChange: c, color: d }) {
    const { classes: h } = M(), [b, P] = (0, E.useLocalStorage)("colorPickerPalette", "set1"), K = m[b], de = Object.keys(m), [ge, he] = (0, l.useState)(d), ve = (0, s.default)(d).rgb().toString(), F = (_) => {
      he(_);
      try {
        c((0, s.default)(_).rgb().toString());
      } catch {
      }
    };
    return l.default.createElement(
      "div",
      { style: { display: "flex", padding: 10 } },
      l.default.createElement(
        "div",
        { style: { width: 200, margin: 5 } },
        l.default.createElement(k.RgbaStringColorPicker, { color: ve, onChange: F })
      ),
      l.default.createElement(
        "div",
        { style: { width: 200, margin: 5 } },
        l.default.createElement(u.Select, { value: b, onChange: (_) => {
          const B = _.target.value;
          P(B);
        } }, de.map((_) => l.default.createElement(u.MenuItem, { value: _, key: _ }, _))),
        l.default.createElement("div", { className: h.swatches }, K.map((_, B) => l.default.createElement("button", { key: `${_}-${B}`, className: h.swatch, style: { background: _ }, onClick: () => F(_) }))),
        l.default.createElement(u.TextField, { helperText: "Manually set color (hex, rgb, or css color name)", value: ge, onChange: (_) => F(_.target.value) })
      )
    );
  }
  e.ColorPicker = g, e.default = e.PopoverPicker;
})(G);
const ut = /* @__PURE__ */ be(G), gt = /* @__PURE__ */ Pe({
  __proto__: null,
  default: ut
}, [G]);
export {
  gt as C,
  G as a,
  ut as b
};
