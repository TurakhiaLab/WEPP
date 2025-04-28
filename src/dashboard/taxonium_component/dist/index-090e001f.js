import { Q as G } from "./index-8ee3baa4.js";
function Q(r, a) {
  for (var e = 0; e < a.length; e++) {
    const n = a[e];
    if (typeof n != "string" && !Array.isArray(n)) {
      for (const t in n)
        if (t !== "default" && !(t in r)) {
          const l = Object.getOwnPropertyDescriptor(n, t);
          l && Object.defineProperty(r, t, l.get ? l : {
            enumerable: !0,
            get: () => n[t]
          });
        }
    }
  }
  return Object.freeze(Object.defineProperty(r, Symbol.toStringTag, { value: "Module" }));
}
var C = { exports: {} }, H = {
  aliceblue: [240, 248, 255],
  antiquewhite: [250, 235, 215],
  aqua: [0, 255, 255],
  aquamarine: [127, 255, 212],
  azure: [240, 255, 255],
  beige: [245, 245, 220],
  bisque: [255, 228, 196],
  black: [0, 0, 0],
  blanchedalmond: [255, 235, 205],
  blue: [0, 0, 255],
  blueviolet: [138, 43, 226],
  brown: [165, 42, 42],
  burlywood: [222, 184, 135],
  cadetblue: [95, 158, 160],
  chartreuse: [127, 255, 0],
  chocolate: [210, 105, 30],
  coral: [255, 127, 80],
  cornflowerblue: [100, 149, 237],
  cornsilk: [255, 248, 220],
  crimson: [220, 20, 60],
  cyan: [0, 255, 255],
  darkblue: [0, 0, 139],
  darkcyan: [0, 139, 139],
  darkgoldenrod: [184, 134, 11],
  darkgray: [169, 169, 169],
  darkgreen: [0, 100, 0],
  darkgrey: [169, 169, 169],
  darkkhaki: [189, 183, 107],
  darkmagenta: [139, 0, 139],
  darkolivegreen: [85, 107, 47],
  darkorange: [255, 140, 0],
  darkorchid: [153, 50, 204],
  darkred: [139, 0, 0],
  darksalmon: [233, 150, 122],
  darkseagreen: [143, 188, 143],
  darkslateblue: [72, 61, 139],
  darkslategray: [47, 79, 79],
  darkslategrey: [47, 79, 79],
  darkturquoise: [0, 206, 209],
  darkviolet: [148, 0, 211],
  deeppink: [255, 20, 147],
  deepskyblue: [0, 191, 255],
  dimgray: [105, 105, 105],
  dimgrey: [105, 105, 105],
  dodgerblue: [30, 144, 255],
  firebrick: [178, 34, 34],
  floralwhite: [255, 250, 240],
  forestgreen: [34, 139, 34],
  fuchsia: [255, 0, 255],
  gainsboro: [220, 220, 220],
  ghostwhite: [248, 248, 255],
  gold: [255, 215, 0],
  goldenrod: [218, 165, 32],
  gray: [128, 128, 128],
  green: [0, 128, 0],
  greenyellow: [173, 255, 47],
  grey: [128, 128, 128],
  honeydew: [240, 255, 240],
  hotpink: [255, 105, 180],
  indianred: [205, 92, 92],
  indigo: [75, 0, 130],
  ivory: [255, 255, 240],
  khaki: [240, 230, 140],
  lavender: [230, 230, 250],
  lavenderblush: [255, 240, 245],
  lawngreen: [124, 252, 0],
  lemonchiffon: [255, 250, 205],
  lightblue: [173, 216, 230],
  lightcoral: [240, 128, 128],
  lightcyan: [224, 255, 255],
  lightgoldenrodyellow: [250, 250, 210],
  lightgray: [211, 211, 211],
  lightgreen: [144, 238, 144],
  lightgrey: [211, 211, 211],
  lightpink: [255, 182, 193],
  lightsalmon: [255, 160, 122],
  lightseagreen: [32, 178, 170],
  lightskyblue: [135, 206, 250],
  lightslategray: [119, 136, 153],
  lightslategrey: [119, 136, 153],
  lightsteelblue: [176, 196, 222],
  lightyellow: [255, 255, 224],
  lime: [0, 255, 0],
  limegreen: [50, 205, 50],
  linen: [250, 240, 230],
  magenta: [255, 0, 255],
  maroon: [128, 0, 0],
  mediumaquamarine: [102, 205, 170],
  mediumblue: [0, 0, 205],
  mediumorchid: [186, 85, 211],
  mediumpurple: [147, 112, 219],
  mediumseagreen: [60, 179, 113],
  mediumslateblue: [123, 104, 238],
  mediumspringgreen: [0, 250, 154],
  mediumturquoise: [72, 209, 204],
  mediumvioletred: [199, 21, 133],
  midnightblue: [25, 25, 112],
  mintcream: [245, 255, 250],
  mistyrose: [255, 228, 225],
  moccasin: [255, 228, 181],
  navajowhite: [255, 222, 173],
  navy: [0, 0, 128],
  oldlace: [253, 245, 230],
  olive: [128, 128, 0],
  olivedrab: [107, 142, 35],
  orange: [255, 165, 0],
  orangered: [255, 69, 0],
  orchid: [218, 112, 214],
  palegoldenrod: [238, 232, 170],
  palegreen: [152, 251, 152],
  paleturquoise: [175, 238, 238],
  palevioletred: [219, 112, 147],
  papayawhip: [255, 239, 213],
  peachpuff: [255, 218, 185],
  peru: [205, 133, 63],
  pink: [255, 192, 203],
  plum: [221, 160, 221],
  powderblue: [176, 224, 230],
  purple: [128, 0, 128],
  rebeccapurple: [102, 51, 153],
  red: [255, 0, 0],
  rosybrown: [188, 143, 143],
  royalblue: [65, 105, 225],
  saddlebrown: [139, 69, 19],
  salmon: [250, 128, 114],
  sandybrown: [244, 164, 96],
  seagreen: [46, 139, 87],
  seashell: [255, 245, 238],
  sienna: [160, 82, 45],
  silver: [192, 192, 192],
  skyblue: [135, 206, 235],
  slateblue: [106, 90, 205],
  slategray: [112, 128, 144],
  slategrey: [112, 128, 144],
  snow: [255, 250, 250],
  springgreen: [0, 255, 127],
  steelblue: [70, 130, 180],
  tan: [210, 180, 140],
  teal: [0, 128, 128],
  thistle: [216, 191, 216],
  tomato: [255, 99, 71],
  turquoise: [64, 224, 208],
  violet: [238, 130, 238],
  wheat: [245, 222, 179],
  white: [255, 255, 255],
  whitesmoke: [245, 245, 245],
  yellow: [255, 255, 0],
  yellowgreen: [154, 205, 50]
}, D = { exports: {} }, V = function(a) {
  return !a || typeof a == "string" ? !1 : a instanceof Array || Array.isArray(a) || a.length >= 0 && (a.splice instanceof Function || Object.getOwnPropertyDescriptor(a, a.length - 1) && a.constructor.name !== "String");
}, W = V, X = Array.prototype.concat, Y = Array.prototype.slice, I = D.exports = function(a) {
  for (var e = [], n = 0, t = a.length; n < t; n++) {
    var l = a[n];
    W(l) ? e = X.call(e, Y.call(l)) : e.push(l);
  }
  return e;
};
I.wrap = function(r) {
  return function() {
    return r(I(arguments));
  };
};
var Z = D.exports, x = H, F = Z, U = Object.hasOwnProperty, _ = /* @__PURE__ */ Object.create(null);
for (var q in x)
  U.call(x, q) && (_[x[q]] = q);
var d = C.exports = {
  to: {},
  get: {}
};
d.get = function(r) {
  var a = r.substring(0, 3).toLowerCase(), e, n;
  switch (a) {
    case "hsl":
      e = d.get.hsl(r), n = "hsl";
      break;
    case "hwb":
      e = d.get.hwb(r), n = "hwb";
      break;
    default:
      e = d.get.rgb(r), n = "rgb";
      break;
  }
  return e ? { model: n, value: e } : null;
};
d.get.rgb = function(r) {
  if (!r)
    return null;
  var a = /^#([a-f0-9]{3,4})$/i, e = /^#([a-f0-9]{6})([a-f0-9]{2})?$/i, n = /^rgba?\(\s*([+-]?\d+)(?=[\s,])\s*(?:,\s*)?([+-]?\d+)(?=[\s,])\s*(?:,\s*)?([+-]?\d+)\s*(?:[,|\/]\s*([+-]?[\d\.]+)(%?)\s*)?\)$/, t = /^rgba?\(\s*([+-]?[\d\.]+)\%\s*,?\s*([+-]?[\d\.]+)\%\s*,?\s*([+-]?[\d\.]+)\%\s*(?:[,|\/]\s*([+-]?[\d\.]+)(%?)\s*)?\)$/, l = /^(\w+)$/, i = [0, 0, 0, 1], o, s, h;
  if (o = r.match(e)) {
    for (h = o[2], o = o[1], s = 0; s < 3; s++) {
      var f = s * 2;
      i[s] = parseInt(o.slice(f, f + 2), 16);
    }
    h && (i[3] = parseInt(h, 16) / 255);
  } else if (o = r.match(a)) {
    for (o = o[1], h = o[3], s = 0; s < 3; s++)
      i[s] = parseInt(o[s] + o[s], 16);
    h && (i[3] = parseInt(h + h, 16) / 255);
  } else if (o = r.match(n)) {
    for (s = 0; s < 3; s++)
      i[s] = parseInt(o[s + 1], 0);
    o[4] && (o[5] ? i[3] = parseFloat(o[4]) * 0.01 : i[3] = parseFloat(o[4]));
  } else if (o = r.match(t)) {
    for (s = 0; s < 3; s++)
      i[s] = Math.round(parseFloat(o[s + 1]) * 2.55);
    o[4] && (o[5] ? i[3] = parseFloat(o[4]) * 0.01 : i[3] = parseFloat(o[4]));
  } else
    return (o = r.match(l)) ? o[1] === "transparent" ? [0, 0, 0, 0] : U.call(x, o[1]) ? (i = x[o[1]], i[3] = 1, i) : null : null;
  for (s = 0; s < 3; s++)
    i[s] = y(i[s], 0, 255);
  return i[3] = y(i[3], 0, 1), i;
};
d.get.hsl = function(r) {
  if (!r)
    return null;
  var a = /^hsla?\(\s*([+-]?(?:\d{0,3}\.)?\d+)(?:deg)?\s*,?\s*([+-]?[\d\.]+)%\s*,?\s*([+-]?[\d\.]+)%\s*(?:[,|\/]\s*([+-]?(?=\.\d|\d)(?:0|[1-9]\d*)?(?:\.\d*)?(?:[eE][+-]?\d+)?)\s*)?\)$/, e = r.match(a);
  if (e) {
    var n = parseFloat(e[4]), t = (parseFloat(e[1]) % 360 + 360) % 360, l = y(parseFloat(e[2]), 0, 100), i = y(parseFloat(e[3]), 0, 100), o = y(isNaN(n) ? 1 : n, 0, 1);
    return [t, l, i, o];
  }
  return null;
};
d.get.hwb = function(r) {
  if (!r)
    return null;
  var a = /^hwb\(\s*([+-]?\d{0,3}(?:\.\d+)?)(?:deg)?\s*,\s*([+-]?[\d\.]+)%\s*,\s*([+-]?[\d\.]+)%\s*(?:,\s*([+-]?(?=\.\d|\d)(?:0|[1-9]\d*)?(?:\.\d*)?(?:[eE][+-]?\d+)?)\s*)?\)$/, e = r.match(a);
  if (e) {
    var n = parseFloat(e[4]), t = (parseFloat(e[1]) % 360 + 360) % 360, l = y(parseFloat(e[2]), 0, 100), i = y(parseFloat(e[3]), 0, 100), o = y(isNaN(n) ? 1 : n, 0, 1);
    return [t, l, i, o];
  }
  return null;
};
d.to.hex = function() {
  var r = F(arguments);
  return "#" + O(r[0]) + O(r[1]) + O(r[2]) + (r[3] < 1 ? O(Math.round(r[3] * 255)) : "");
};
d.to.rgb = function() {
  var r = F(arguments);
  return r.length < 4 || r[3] === 1 ? "rgb(" + Math.round(r[0]) + ", " + Math.round(r[1]) + ", " + Math.round(r[2]) + ")" : "rgba(" + Math.round(r[0]) + ", " + Math.round(r[1]) + ", " + Math.round(r[2]) + ", " + r[3] + ")";
};
d.to.rgb.percent = function() {
  var r = F(arguments), a = Math.round(r[0] / 255 * 100), e = Math.round(r[1] / 255 * 100), n = Math.round(r[2] / 255 * 100);
  return r.length < 4 || r[3] === 1 ? "rgb(" + a + "%, " + e + "%, " + n + "%)" : "rgba(" + a + "%, " + e + "%, " + n + "%, " + r[3] + ")";
};
d.to.hsl = function() {
  var r = F(arguments);
  return r.length < 4 || r[3] === 1 ? "hsl(" + r[0] + ", " + r[1] + "%, " + r[2] + "%)" : "hsla(" + r[0] + ", " + r[1] + "%, " + r[2] + "%, " + r[3] + ")";
};
d.to.hwb = function() {
  var r = F(arguments), a = "";
  return r.length >= 4 && r[3] !== 1 && (a = ", " + r[3]), "hwb(" + r[0] + ", " + r[1] + "%, " + r[2] + "%" + a + ")";
};
d.to.keyword = function(r) {
  return _[r.slice(0, 3)];
};
function y(r, a, e) {
  return Math.min(Math.max(a, r), e);
}
function O(r) {
  var a = Math.round(r).toString(16).toUpperCase();
  return a.length < 2 ? "0" + a : a;
}
var j = C.exports, K = { exports: {} }, rr = {
  aliceblue: [240, 248, 255],
  antiquewhite: [250, 235, 215],
  aqua: [0, 255, 255],
  aquamarine: [127, 255, 212],
  azure: [240, 255, 255],
  beige: [245, 245, 220],
  bisque: [255, 228, 196],
  black: [0, 0, 0],
  blanchedalmond: [255, 235, 205],
  blue: [0, 0, 255],
  blueviolet: [138, 43, 226],
  brown: [165, 42, 42],
  burlywood: [222, 184, 135],
  cadetblue: [95, 158, 160],
  chartreuse: [127, 255, 0],
  chocolate: [210, 105, 30],
  coral: [255, 127, 80],
  cornflowerblue: [100, 149, 237],
  cornsilk: [255, 248, 220],
  crimson: [220, 20, 60],
  cyan: [0, 255, 255],
  darkblue: [0, 0, 139],
  darkcyan: [0, 139, 139],
  darkgoldenrod: [184, 134, 11],
  darkgray: [169, 169, 169],
  darkgreen: [0, 100, 0],
  darkgrey: [169, 169, 169],
  darkkhaki: [189, 183, 107],
  darkmagenta: [139, 0, 139],
  darkolivegreen: [85, 107, 47],
  darkorange: [255, 140, 0],
  darkorchid: [153, 50, 204],
  darkred: [139, 0, 0],
  darksalmon: [233, 150, 122],
  darkseagreen: [143, 188, 143],
  darkslateblue: [72, 61, 139],
  darkslategray: [47, 79, 79],
  darkslategrey: [47, 79, 79],
  darkturquoise: [0, 206, 209],
  darkviolet: [148, 0, 211],
  deeppink: [255, 20, 147],
  deepskyblue: [0, 191, 255],
  dimgray: [105, 105, 105],
  dimgrey: [105, 105, 105],
  dodgerblue: [30, 144, 255],
  firebrick: [178, 34, 34],
  floralwhite: [255, 250, 240],
  forestgreen: [34, 139, 34],
  fuchsia: [255, 0, 255],
  gainsboro: [220, 220, 220],
  ghostwhite: [248, 248, 255],
  gold: [255, 215, 0],
  goldenrod: [218, 165, 32],
  gray: [128, 128, 128],
  green: [0, 128, 0],
  greenyellow: [173, 255, 47],
  grey: [128, 128, 128],
  honeydew: [240, 255, 240],
  hotpink: [255, 105, 180],
  indianred: [205, 92, 92],
  indigo: [75, 0, 130],
  ivory: [255, 255, 240],
  khaki: [240, 230, 140],
  lavender: [230, 230, 250],
  lavenderblush: [255, 240, 245],
  lawngreen: [124, 252, 0],
  lemonchiffon: [255, 250, 205],
  lightblue: [173, 216, 230],
  lightcoral: [240, 128, 128],
  lightcyan: [224, 255, 255],
  lightgoldenrodyellow: [250, 250, 210],
  lightgray: [211, 211, 211],
  lightgreen: [144, 238, 144],
  lightgrey: [211, 211, 211],
  lightpink: [255, 182, 193],
  lightsalmon: [255, 160, 122],
  lightseagreen: [32, 178, 170],
  lightskyblue: [135, 206, 250],
  lightslategray: [119, 136, 153],
  lightslategrey: [119, 136, 153],
  lightsteelblue: [176, 196, 222],
  lightyellow: [255, 255, 224],
  lime: [0, 255, 0],
  limegreen: [50, 205, 50],
  linen: [250, 240, 230],
  magenta: [255, 0, 255],
  maroon: [128, 0, 0],
  mediumaquamarine: [102, 205, 170],
  mediumblue: [0, 0, 205],
  mediumorchid: [186, 85, 211],
  mediumpurple: [147, 112, 219],
  mediumseagreen: [60, 179, 113],
  mediumslateblue: [123, 104, 238],
  mediumspringgreen: [0, 250, 154],
  mediumturquoise: [72, 209, 204],
  mediumvioletred: [199, 21, 133],
  midnightblue: [25, 25, 112],
  mintcream: [245, 255, 250],
  mistyrose: [255, 228, 225],
  moccasin: [255, 228, 181],
  navajowhite: [255, 222, 173],
  navy: [0, 0, 128],
  oldlace: [253, 245, 230],
  olive: [128, 128, 0],
  olivedrab: [107, 142, 35],
  orange: [255, 165, 0],
  orangered: [255, 69, 0],
  orchid: [218, 112, 214],
  palegoldenrod: [238, 232, 170],
  palegreen: [152, 251, 152],
  paleturquoise: [175, 238, 238],
  palevioletred: [219, 112, 147],
  papayawhip: [255, 239, 213],
  peachpuff: [255, 218, 185],
  peru: [205, 133, 63],
  pink: [255, 192, 203],
  plum: [221, 160, 221],
  powderblue: [176, 224, 230],
  purple: [128, 0, 128],
  rebeccapurple: [102, 51, 153],
  red: [255, 0, 0],
  rosybrown: [188, 143, 143],
  royalblue: [65, 105, 225],
  saddlebrown: [139, 69, 19],
  salmon: [250, 128, 114],
  sandybrown: [244, 164, 96],
  seagreen: [46, 139, 87],
  seashell: [255, 245, 238],
  sienna: [160, 82, 45],
  silver: [192, 192, 192],
  skyblue: [135, 206, 235],
  slateblue: [106, 90, 205],
  slategray: [112, 128, 144],
  slategrey: [112, 128, 144],
  snow: [255, 250, 250],
  springgreen: [0, 255, 127],
  steelblue: [70, 130, 180],
  tan: [210, 180, 140],
  teal: [0, 128, 128],
  thistle: [216, 191, 216],
  tomato: [255, 99, 71],
  turquoise: [64, 224, 208],
  violet: [238, 130, 238],
  wheat: [245, 222, 179],
  white: [255, 255, 255],
  whitesmoke: [245, 245, 245],
  yellow: [255, 255, 0],
  yellowgreen: [154, 205, 50]
}, w = rr, R = {};
for (var E in w)
  w.hasOwnProperty(E) && (R[w[E]] = E);
var u = K.exports = {
  rgb: { channels: 3, labels: "rgb" },
  hsl: { channels: 3, labels: "hsl" },
  hsv: { channels: 3, labels: "hsv" },
  hwb: { channels: 3, labels: "hwb" },
  cmyk: { channels: 4, labels: "cmyk" },
  xyz: { channels: 3, labels: "xyz" },
  lab: { channels: 3, labels: "lab" },
  lch: { channels: 3, labels: "lch" },
  hex: { channels: 1, labels: ["hex"] },
  keyword: { channels: 1, labels: ["keyword"] },
  ansi16: { channels: 1, labels: ["ansi16"] },
  ansi256: { channels: 1, labels: ["ansi256"] },
  hcg: { channels: 3, labels: ["h", "c", "g"] },
  apple: { channels: 3, labels: ["r16", "g16", "b16"] },
  gray: { channels: 1, labels: ["gray"] }
};
for (var b in u)
  if (u.hasOwnProperty(b)) {
    if (!("channels" in u[b]))
      throw new Error("missing channels property: " + b);
    if (!("labels" in u[b]))
      throw new Error("missing channel labels property: " + b);
    if (u[b].labels.length !== u[b].channels)
      throw new Error("channel and label counts mismatch: " + b);
    var ar = u[b].channels, er = u[b].labels;
    delete u[b].channels, delete u[b].labels, Object.defineProperty(u[b], "channels", { value: ar }), Object.defineProperty(u[b], "labels", { value: er });
  }
u.rgb.hsl = function(r) {
  var a = r[0] / 255, e = r[1] / 255, n = r[2] / 255, t = Math.min(a, e, n), l = Math.max(a, e, n), i = l - t, o, s, h;
  return l === t ? o = 0 : a === l ? o = (e - n) / i : e === l ? o = 2 + (n - a) / i : n === l && (o = 4 + (a - e) / i), o = Math.min(o * 60, 360), o < 0 && (o += 360), h = (t + l) / 2, l === t ? s = 0 : h <= 0.5 ? s = i / (l + t) : s = i / (2 - l - t), [o, s * 100, h * 100];
};
u.rgb.hsv = function(r) {
  var a, e, n, t, l, i = r[0] / 255, o = r[1] / 255, s = r[2] / 255, h = Math.max(i, o, s), f = h - Math.min(i, o, s), m = function(B) {
    return (h - B) / 6 / f + 1 / 2;
  };
  return f === 0 ? t = l = 0 : (l = f / h, a = m(i), e = m(o), n = m(s), i === h ? t = n - e : o === h ? t = 1 / 3 + a - n : s === h && (t = 2 / 3 + e - a), t < 0 ? t += 1 : t > 1 && (t -= 1)), [
    t * 360,
    l * 100,
    h * 100
  ];
};
u.rgb.hwb = function(r) {
  var a = r[0], e = r[1], n = r[2], t = u.rgb.hsl(r)[0], l = 1 / 255 * Math.min(a, Math.min(e, n));
  return n = 1 - 1 / 255 * Math.max(a, Math.max(e, n)), [t, l * 100, n * 100];
};
u.rgb.cmyk = function(r) {
  var a = r[0] / 255, e = r[1] / 255, n = r[2] / 255, t, l, i, o;
  return o = Math.min(1 - a, 1 - e, 1 - n), t = (1 - a - o) / (1 - o) || 0, l = (1 - e - o) / (1 - o) || 0, i = (1 - n - o) / (1 - o) || 0, [t * 100, l * 100, i * 100, o * 100];
};
function nr(r, a) {
  return Math.pow(r[0] - a[0], 2) + Math.pow(r[1] - a[1], 2) + Math.pow(r[2] - a[2], 2);
}
u.rgb.keyword = function(r) {
  var a = R[r];
  if (a)
    return a;
  var e = 1 / 0, n;
  for (var t in w)
    if (w.hasOwnProperty(t)) {
      var l = w[t], i = nr(r, l);
      i < e && (e = i, n = t);
    }
  return n;
};
u.keyword.rgb = function(r) {
  return w[r];
};
u.rgb.xyz = function(r) {
  var a = r[0] / 255, e = r[1] / 255, n = r[2] / 255;
  a = a > 0.04045 ? Math.pow((a + 0.055) / 1.055, 2.4) : a / 12.92, e = e > 0.04045 ? Math.pow((e + 0.055) / 1.055, 2.4) : e / 12.92, n = n > 0.04045 ? Math.pow((n + 0.055) / 1.055, 2.4) : n / 12.92;
  var t = a * 0.4124 + e * 0.3576 + n * 0.1805, l = a * 0.2126 + e * 0.7152 + n * 0.0722, i = a * 0.0193 + e * 0.1192 + n * 0.9505;
  return [t * 100, l * 100, i * 100];
};
u.rgb.lab = function(r) {
  var a = u.rgb.xyz(r), e = a[0], n = a[1], t = a[2], l, i, o;
  return e /= 95.047, n /= 100, t /= 108.883, e = e > 8856e-6 ? Math.pow(e, 1 / 3) : 7.787 * e + 16 / 116, n = n > 8856e-6 ? Math.pow(n, 1 / 3) : 7.787 * n + 16 / 116, t = t > 8856e-6 ? Math.pow(t, 1 / 3) : 7.787 * t + 16 / 116, l = 116 * n - 16, i = 500 * (e - n), o = 200 * (n - t), [l, i, o];
};
u.hsl.rgb = function(r) {
  var a = r[0] / 360, e = r[1] / 100, n = r[2] / 100, t, l, i, o, s;
  if (e === 0)
    return s = n * 255, [s, s, s];
  n < 0.5 ? l = n * (1 + e) : l = n + e - n * e, t = 2 * n - l, o = [0, 0, 0];
  for (var h = 0; h < 3; h++)
    i = a + 1 / 3 * -(h - 1), i < 0 && i++, i > 1 && i--, 6 * i < 1 ? s = t + (l - t) * 6 * i : 2 * i < 1 ? s = l : 3 * i < 2 ? s = t + (l - t) * (2 / 3 - i) * 6 : s = t, o[h] = s * 255;
  return o;
};
u.hsl.hsv = function(r) {
  var a = r[0], e = r[1] / 100, n = r[2] / 100, t = e, l = Math.max(n, 0.01), i, o;
  return n *= 2, e *= n <= 1 ? n : 2 - n, t *= l <= 1 ? l : 2 - l, o = (n + e) / 2, i = n === 0 ? 2 * t / (l + t) : 2 * e / (n + e), [a, i * 100, o * 100];
};
u.hsv.rgb = function(r) {
  var a = r[0] / 60, e = r[1] / 100, n = r[2] / 100, t = Math.floor(a) % 6, l = a - Math.floor(a), i = 255 * n * (1 - e), o = 255 * n * (1 - e * l), s = 255 * n * (1 - e * (1 - l));
  switch (n *= 255, t) {
    case 0:
      return [n, s, i];
    case 1:
      return [o, n, i];
    case 2:
      return [i, n, s];
    case 3:
      return [i, o, n];
    case 4:
      return [s, i, n];
    case 5:
      return [n, i, o];
  }
};
u.hsv.hsl = function(r) {
  var a = r[0], e = r[1] / 100, n = r[2] / 100, t = Math.max(n, 0.01), l, i, o;
  return o = (2 - e) * n, l = (2 - e) * t, i = e * t, i /= l <= 1 ? l : 2 - l, i = i || 0, o /= 2, [a, i * 100, o * 100];
};
u.hwb.rgb = function(r) {
  var a = r[0] / 360, e = r[1] / 100, n = r[2] / 100, t = e + n, l, i, o, s;
  t > 1 && (e /= t, n /= t), l = Math.floor(6 * a), i = 1 - n, o = 6 * a - l, l & 1 && (o = 1 - o), s = e + o * (i - e);
  var h, f, m;
  switch (l) {
    default:
    case 6:
    case 0:
      h = i, f = s, m = e;
      break;
    case 1:
      h = s, f = i, m = e;
      break;
    case 2:
      h = e, f = i, m = s;
      break;
    case 3:
      h = e, f = s, m = i;
      break;
    case 4:
      h = s, f = e, m = i;
      break;
    case 5:
      h = i, f = e, m = s;
      break;
  }
  return [h * 255, f * 255, m * 255];
};
u.cmyk.rgb = function(r) {
  var a = r[0] / 100, e = r[1] / 100, n = r[2] / 100, t = r[3] / 100, l, i, o;
  return l = 1 - Math.min(1, a * (1 - t) + t), i = 1 - Math.min(1, e * (1 - t) + t), o = 1 - Math.min(1, n * (1 - t) + t), [l * 255, i * 255, o * 255];
};
u.xyz.rgb = function(r) {
  var a = r[0] / 100, e = r[1] / 100, n = r[2] / 100, t, l, i;
  return t = a * 3.2406 + e * -1.5372 + n * -0.4986, l = a * -0.9689 + e * 1.8758 + n * 0.0415, i = a * 0.0557 + e * -0.204 + n * 1.057, t = t > 31308e-7 ? 1.055 * Math.pow(t, 1 / 2.4) - 0.055 : t * 12.92, l = l > 31308e-7 ? 1.055 * Math.pow(l, 1 / 2.4) - 0.055 : l * 12.92, i = i > 31308e-7 ? 1.055 * Math.pow(i, 1 / 2.4) - 0.055 : i * 12.92, t = Math.min(Math.max(0, t), 1), l = Math.min(Math.max(0, l), 1), i = Math.min(Math.max(0, i), 1), [t * 255, l * 255, i * 255];
};
u.xyz.lab = function(r) {
  var a = r[0], e = r[1], n = r[2], t, l, i;
  return a /= 95.047, e /= 100, n /= 108.883, a = a > 8856e-6 ? Math.pow(a, 1 / 3) : 7.787 * a + 16 / 116, e = e > 8856e-6 ? Math.pow(e, 1 / 3) : 7.787 * e + 16 / 116, n = n > 8856e-6 ? Math.pow(n, 1 / 3) : 7.787 * n + 16 / 116, t = 116 * e - 16, l = 500 * (a - e), i = 200 * (e - n), [t, l, i];
};
u.lab.xyz = function(r) {
  var a = r[0], e = r[1], n = r[2], t, l, i;
  l = (a + 16) / 116, t = e / 500 + l, i = l - n / 200;
  var o = Math.pow(l, 3), s = Math.pow(t, 3), h = Math.pow(i, 3);
  return l = o > 8856e-6 ? o : (l - 16 / 116) / 7.787, t = s > 8856e-6 ? s : (t - 16 / 116) / 7.787, i = h > 8856e-6 ? h : (i - 16 / 116) / 7.787, t *= 95.047, l *= 100, i *= 108.883, [t, l, i];
};
u.lab.lch = function(r) {
  var a = r[0], e = r[1], n = r[2], t, l, i;
  return t = Math.atan2(n, e), l = t * 360 / 2 / Math.PI, l < 0 && (l += 360), i = Math.sqrt(e * e + n * n), [a, i, l];
};
u.lch.lab = function(r) {
  var a = r[0], e = r[1], n = r[2], t, l, i;
  return i = n / 360 * 2 * Math.PI, t = e * Math.cos(i), l = e * Math.sin(i), [a, t, l];
};
u.rgb.ansi16 = function(r) {
  var a = r[0], e = r[1], n = r[2], t = 1 in arguments ? arguments[1] : u.rgb.hsv(r)[2];
  if (t = Math.round(t / 50), t === 0)
    return 30;
  var l = 30 + (Math.round(n / 255) << 2 | Math.round(e / 255) << 1 | Math.round(a / 255));
  return t === 2 && (l += 60), l;
};
u.hsv.ansi16 = function(r) {
  return u.rgb.ansi16(u.hsv.rgb(r), r[2]);
};
u.rgb.ansi256 = function(r) {
  var a = r[0], e = r[1], n = r[2];
  if (a === e && e === n)
    return a < 8 ? 16 : a > 248 ? 231 : Math.round((a - 8) / 247 * 24) + 232;
  var t = 16 + 36 * Math.round(a / 255 * 5) + 6 * Math.round(e / 255 * 5) + Math.round(n / 255 * 5);
  return t;
};
u.ansi16.rgb = function(r) {
  var a = r % 10;
  if (a === 0 || a === 7)
    return r > 50 && (a += 3.5), a = a / 10.5 * 255, [a, a, a];
  var e = (~~(r > 50) + 1) * 0.5, n = (a & 1) * e * 255, t = (a >> 1 & 1) * e * 255, l = (a >> 2 & 1) * e * 255;
  return [n, t, l];
};
u.ansi256.rgb = function(r) {
  if (r >= 232) {
    var a = (r - 232) * 10 + 8;
    return [a, a, a];
  }
  r -= 16;
  var e, n = Math.floor(r / 36) / 5 * 255, t = Math.floor((e = r % 36) / 6) / 5 * 255, l = e % 6 / 5 * 255;
  return [n, t, l];
};
u.rgb.hex = function(r) {
  var a = ((Math.round(r[0]) & 255) << 16) + ((Math.round(r[1]) & 255) << 8) + (Math.round(r[2]) & 255), e = a.toString(16).toUpperCase();
  return "000000".substring(e.length) + e;
};
u.hex.rgb = function(r) {
  var a = r.toString(16).match(/[a-f0-9]{6}|[a-f0-9]{3}/i);
  if (!a)
    return [0, 0, 0];
  var e = a[0];
  a[0].length === 3 && (e = e.split("").map(function(o) {
    return o + o;
  }).join(""));
  var n = parseInt(e, 16), t = n >> 16 & 255, l = n >> 8 & 255, i = n & 255;
  return [t, l, i];
};
u.rgb.hcg = function(r) {
  var a = r[0] / 255, e = r[1] / 255, n = r[2] / 255, t = Math.max(Math.max(a, e), n), l = Math.min(Math.min(a, e), n), i = t - l, o, s;
  return i < 1 ? o = l / (1 - i) : o = 0, i <= 0 ? s = 0 : t === a ? s = (e - n) / i % 6 : t === e ? s = 2 + (n - a) / i : s = 4 + (a - e) / i + 4, s /= 6, s %= 1, [s * 360, i * 100, o * 100];
};
u.hsl.hcg = function(r) {
  var a = r[1] / 100, e = r[2] / 100, n = 1, t = 0;
  return e < 0.5 ? n = 2 * a * e : n = 2 * a * (1 - e), n < 1 && (t = (e - 0.5 * n) / (1 - n)), [r[0], n * 100, t * 100];
};
u.hsv.hcg = function(r) {
  var a = r[1] / 100, e = r[2] / 100, n = a * e, t = 0;
  return n < 1 && (t = (e - n) / (1 - n)), [r[0], n * 100, t * 100];
};
u.hcg.rgb = function(r) {
  var a = r[0] / 360, e = r[1] / 100, n = r[2] / 100;
  if (e === 0)
    return [n * 255, n * 255, n * 255];
  var t = [0, 0, 0], l = a % 1 * 6, i = l % 1, o = 1 - i, s = 0;
  switch (Math.floor(l)) {
    case 0:
      t[0] = 1, t[1] = i, t[2] = 0;
      break;
    case 1:
      t[0] = o, t[1] = 1, t[2] = 0;
      break;
    case 2:
      t[0] = 0, t[1] = 1, t[2] = i;
      break;
    case 3:
      t[0] = 0, t[1] = o, t[2] = 1;
      break;
    case 4:
      t[0] = i, t[1] = 0, t[2] = 1;
      break;
    default:
      t[0] = 1, t[1] = 0, t[2] = o;
  }
  return s = (1 - e) * n, [
    (e * t[0] + s) * 255,
    (e * t[1] + s) * 255,
    (e * t[2] + s) * 255
  ];
};
u.hcg.hsv = function(r) {
  var a = r[1] / 100, e = r[2] / 100, n = a + e * (1 - a), t = 0;
  return n > 0 && (t = a / n), [r[0], t * 100, n * 100];
};
u.hcg.hsl = function(r) {
  var a = r[1] / 100, e = r[2] / 100, n = e * (1 - a) + 0.5 * a, t = 0;
  return n > 0 && n < 0.5 ? t = a / (2 * n) : n >= 0.5 && n < 1 && (t = a / (2 * (1 - n))), [r[0], t * 100, n * 100];
};
u.hcg.hwb = function(r) {
  var a = r[1] / 100, e = r[2] / 100, n = a + e * (1 - a);
  return [r[0], (n - a) * 100, (1 - n) * 100];
};
u.hwb.hcg = function(r) {
  var a = r[1] / 100, e = r[2] / 100, n = 1 - e, t = n - a, l = 0;
  return t < 1 && (l = (n - t) / (1 - t)), [r[0], t * 100, l * 100];
};
u.apple.rgb = function(r) {
  return [r[0] / 65535 * 255, r[1] / 65535 * 255, r[2] / 65535 * 255];
};
u.rgb.apple = function(r) {
  return [r[0] / 255 * 65535, r[1] / 255 * 65535, r[2] / 255 * 65535];
};
u.gray.rgb = function(r) {
  return [r[0] / 100 * 255, r[0] / 100 * 255, r[0] / 100 * 255];
};
u.gray.hsl = u.gray.hsv = function(r) {
  return [0, 0, r[0]];
};
u.gray.hwb = function(r) {
  return [0, 100, r[0]];
};
u.gray.cmyk = function(r) {
  return [0, 0, 0, r[0]];
};
u.gray.lab = function(r) {
  return [r[0], 0, 0];
};
u.gray.hex = function(r) {
  var a = Math.round(r[0] / 100 * 255) & 255, e = (a << 16) + (a << 8) + a, n = e.toString(16).toUpperCase();
  return "000000".substring(n.length) + n;
};
u.rgb.gray = function(r) {
  var a = (r[0] + r[1] + r[2]) / 3;
  return [a / 255 * 100];
};
var T = K.exports, z = T;
function tr() {
  for (var r = {}, a = Object.keys(z), e = a.length, n = 0; n < e; n++)
    r[a[n]] = {
      // http://jsperf.com/1-vs-infinity
      // micro-opt, but this is simple.
      distance: -1,
      parent: null
    };
  return r;
}
function ir(r) {
  var a = tr(), e = [r];
  for (a[r].distance = 0; e.length; )
    for (var n = e.pop(), t = Object.keys(z[n]), l = t.length, i = 0; i < l; i++) {
      var o = t[i], s = a[o];
      s.distance === -1 && (s.distance = a[n].distance + 1, s.parent = n, e.unshift(o));
    }
  return a;
}
function lr(r, a) {
  return function(e) {
    return a(r(e));
  };
}
function or(r, a) {
  for (var e = [a[r].parent, r], n = z[a[r].parent][r], t = a[r].parent; a[t].parent; )
    e.unshift(a[t].parent), n = lr(z[a[t].parent][t], n), t = a[t].parent;
  return n.conversion = e, n;
}
var sr = function(r) {
  for (var a = ir(r), e = {}, n = Object.keys(a), t = n.length, l = 0; l < t; l++) {
    var i = n[l], o = a[i];
    o.parent !== null && (e[i] = or(i, a));
  }
  return e;
}, S = T, ur = sr, k = {}, hr = Object.keys(S);
function vr(r) {
  var a = function(e) {
    return e == null ? e : (arguments.length > 1 && (e = Array.prototype.slice.call(arguments)), r(e));
  };
  return "conversion" in r && (a.conversion = r.conversion), a;
}
function cr(r) {
  var a = function(e) {
    if (e == null)
      return e;
    arguments.length > 1 && (e = Array.prototype.slice.call(arguments));
    var n = r(e);
    if (typeof n == "object")
      for (var t = n.length, l = 0; l < t; l++)
        n[l] = Math.round(n[l]);
    return n;
  };
  return "conversion" in r && (a.conversion = r.conversion), a;
}
hr.forEach(function(r) {
  k[r] = {}, Object.defineProperty(k[r], "channels", { value: S[r].channels }), Object.defineProperty(k[r], "labels", { value: S[r].labels });
  var a = ur(r), e = Object.keys(a);
  e.forEach(function(n) {
    var t = a[n];
    k[r][n] = cr(t), k[r][n].raw = vr(t);
  });
});
var fr = k, M = j, p = fr, N = [].slice, J = [
  // to be honest, I don't really feel like keyword belongs in color convert, but eh.
  "keyword",
  // gray conflicts with some method names, and has its own method defined.
  "gray",
  // shouldn't really be in color-convert either...
  "hex"
], $ = {};
Object.keys(p).forEach(function(r) {
  $[N.call(p[r].labels).sort().join("")] = r;
});
var A = {};
function g(r, a) {
  if (!(this instanceof g))
    return new g(r, a);
  if (a && a in J && (a = null), a && !(a in p))
    throw new Error("Unknown model: " + a);
  var e, n;
  if (r == null)
    this.model = "rgb", this.color = [0, 0, 0], this.valpha = 1;
  else if (r instanceof g)
    this.model = r.model, this.color = r.color.slice(), this.valpha = r.valpha;
  else if (typeof r == "string") {
    var t = M.get(r);
    if (t === null)
      throw new Error("Unable to parse color from string: " + r);
    this.model = t.model, n = p[this.model].channels, this.color = t.value.slice(0, n), this.valpha = typeof t.value[n] == "number" ? t.value[n] : 1;
  } else if (r.length) {
    this.model = a || "rgb", n = p[this.model].channels;
    var l = N.call(r, 0, n);
    this.color = P(l, n), this.valpha = typeof r[n] == "number" ? r[n] : 1;
  } else if (typeof r == "number")
    r &= 16777215, this.model = "rgb", this.color = [
      r >> 16 & 255,
      r >> 8 & 255,
      r & 255
    ], this.valpha = 1;
  else {
    this.valpha = 1;
    var i = Object.keys(r);
    "alpha" in r && (i.splice(i.indexOf("alpha"), 1), this.valpha = typeof r.alpha == "number" ? r.alpha : 0);
    var o = i.sort().join("");
    if (!(o in $))
      throw new Error("Unable to parse color from object: " + JSON.stringify(r));
    this.model = $[o];
    var s = p[this.model].labels, h = [];
    for (e = 0; e < s.length; e++)
      h.push(r[s[e]]);
    this.color = P(h);
  }
  if (A[this.model])
    for (n = p[this.model].channels, e = 0; e < n; e++) {
      var f = A[this.model][e];
      f && (this.color[e] = f(this.color[e]));
    }
  this.valpha = Math.max(0, Math.min(1, this.valpha)), Object.freeze && Object.freeze(this);
}
g.prototype = {
  toString: function() {
    return this.string();
  },
  toJSON: function() {
    return this[this.model]();
  },
  string: function(r) {
    var a = this.model in M.to ? this : this.rgb();
    a = a.round(typeof r == "number" ? r : 1);
    var e = a.valpha === 1 ? a.color : a.color.concat(this.valpha);
    return M.to[a.model](e);
  },
  percentString: function(r) {
    var a = this.rgb().round(typeof r == "number" ? r : 1), e = a.valpha === 1 ? a.color : a.color.concat(this.valpha);
    return M.to.rgb.percent(e);
  },
  array: function() {
    return this.valpha === 1 ? this.color.slice() : this.color.concat(this.valpha);
  },
  object: function() {
    for (var r = {}, a = p[this.model].channels, e = p[this.model].labels, n = 0; n < a; n++)
      r[e[n]] = this.color[n];
    return this.valpha !== 1 && (r.alpha = this.valpha), r;
  },
  unitArray: function() {
    var r = this.rgb().color;
    return r[0] /= 255, r[1] /= 255, r[2] /= 255, this.valpha !== 1 && r.push(this.valpha), r;
  },
  unitObject: function() {
    var r = this.rgb().object();
    return r.r /= 255, r.g /= 255, r.b /= 255, this.valpha !== 1 && (r.alpha = this.valpha), r;
  },
  round: function(r) {
    return r = Math.max(r || 0, 0), new g(this.color.map(br(r)).concat(this.valpha), this.model);
  },
  alpha: function(r) {
    return arguments.length ? new g(this.color.concat(Math.max(0, Math.min(1, r))), this.model) : this.valpha;
  },
  // rgb
  red: v("rgb", 0, c(255)),
  green: v("rgb", 1, c(255)),
  blue: v("rgb", 2, c(255)),
  hue: v(["hsl", "hsv", "hsl", "hwb", "hcg"], 0, function(r) {
    return (r % 360 + 360) % 360;
  }),
  // eslint-disable-line brace-style
  saturationl: v("hsl", 1, c(100)),
  lightness: v("hsl", 2, c(100)),
  saturationv: v("hsv", 1, c(100)),
  value: v("hsv", 2, c(100)),
  chroma: v("hcg", 1, c(100)),
  gray: v("hcg", 2, c(100)),
  white: v("hwb", 1, c(100)),
  wblack: v("hwb", 2, c(100)),
  cyan: v("cmyk", 0, c(100)),
  magenta: v("cmyk", 1, c(100)),
  yellow: v("cmyk", 2, c(100)),
  black: v("cmyk", 3, c(100)),
  x: v("xyz", 0, c(100)),
  y: v("xyz", 1, c(100)),
  z: v("xyz", 2, c(100)),
  l: v("lab", 0, c(100)),
  a: v("lab", 1),
  b: v("lab", 2),
  keyword: function(r) {
    return arguments.length ? new g(r) : p[this.model].keyword(this.color);
  },
  hex: function(r) {
    return arguments.length ? new g(r) : M.to.hex(this.rgb().round().color);
  },
  rgbNumber: function() {
    var r = this.rgb().color;
    return (r[0] & 255) << 16 | (r[1] & 255) << 8 | r[2] & 255;
  },
  luminosity: function() {
    for (var r = this.rgb().color, a = [], e = 0; e < r.length; e++) {
      var n = r[e] / 255;
      a[e] = n <= 0.03928 ? n / 12.92 : Math.pow((n + 0.055) / 1.055, 2.4);
    }
    return 0.2126 * a[0] + 0.7152 * a[1] + 0.0722 * a[2];
  },
  contrast: function(r) {
    var a = this.luminosity(), e = r.luminosity();
    return a > e ? (a + 0.05) / (e + 0.05) : (e + 0.05) / (a + 0.05);
  },
  level: function(r) {
    var a = this.contrast(r);
    return a >= 7.1 ? "AAA" : a >= 4.5 ? "AA" : "";
  },
  isDark: function() {
    var r = this.rgb().color, a = (r[0] * 299 + r[1] * 587 + r[2] * 114) / 1e3;
    return a < 128;
  },
  isLight: function() {
    return !this.isDark();
  },
  negate: function() {
    for (var r = this.rgb(), a = 0; a < 3; a++)
      r.color[a] = 255 - r.color[a];
    return r;
  },
  lighten: function(r) {
    var a = this.hsl();
    return a.color[2] += a.color[2] * r, a;
  },
  darken: function(r) {
    var a = this.hsl();
    return a.color[2] -= a.color[2] * r, a;
  },
  saturate: function(r) {
    var a = this.hsl();
    return a.color[1] += a.color[1] * r, a;
  },
  desaturate: function(r) {
    var a = this.hsl();
    return a.color[1] -= a.color[1] * r, a;
  },
  whiten: function(r) {
    var a = this.hwb();
    return a.color[1] += a.color[1] * r, a;
  },
  blacken: function(r) {
    var a = this.hwb();
    return a.color[2] += a.color[2] * r, a;
  },
  grayscale: function() {
    var r = this.rgb().color, a = r[0] * 0.3 + r[1] * 0.59 + r[2] * 0.11;
    return g.rgb(a, a, a);
  },
  fade: function(r) {
    return this.alpha(this.valpha - this.valpha * r);
  },
  opaquer: function(r) {
    return this.alpha(this.valpha + this.valpha * r);
  },
  rotate: function(r) {
    var a = this.hsl(), e = a.color[0];
    return e = (e + r) % 360, e = e < 0 ? 360 + e : e, a.color[0] = e, a;
  },
  mix: function(r, a) {
    if (!r || !r.rgb)
      throw new Error('Argument to "mix" was not a Color instance, but rather an instance of ' + typeof r);
    var e = r.rgb(), n = this.rgb(), t = a === void 0 ? 0.5 : a, l = 2 * t - 1, i = e.alpha() - n.alpha(), o = ((l * i === -1 ? l : (l + i) / (1 + l * i)) + 1) / 2, s = 1 - o;
    return g.rgb(
      o * e.red() + s * n.red(),
      o * e.green() + s * n.green(),
      o * e.blue() + s * n.blue(),
      e.alpha() * t + n.alpha() * (1 - t)
    );
  }
};
Object.keys(p).forEach(function(r) {
  if (J.indexOf(r) === -1) {
    var a = p[r].channels;
    g.prototype[r] = function() {
      if (this.model === r)
        return new g(this);
      if (arguments.length)
        return new g(arguments, r);
      var e = typeof arguments[a] == "number" ? a : this.valpha;
      return new g(pr(p[this.model][r].raw(this.color)).concat(e), r);
    }, g[r] = function(e) {
      return typeof e == "number" && (e = P(N.call(arguments), a)), new g(e, r);
    };
  }
});
function gr(r, a) {
  return Number(r.toFixed(a));
}
function br(r) {
  return function(a) {
    return gr(a, r);
  };
}
function v(r, a, e) {
  return r = Array.isArray(r) ? r : [r], r.forEach(function(n) {
    (A[n] || (A[n] = []))[a] = e;
  }), r = r[0], function(n) {
    var t;
    return arguments.length ? (e && (n = e(n)), t = this[r](), t.color[a] = n, t) : (t = this[r]().color[a], e && (t = e(t)), t);
  };
}
function c(r) {
  return function(a) {
    return Math.max(0, Math.min(r, a));
  };
}
function pr(r) {
  return Array.isArray(r) ? r : [r];
}
function P(r, a) {
  for (var e = 0; e < a; e++)
    typeof r[e] != "number" && (r[e] = 0);
  return r;
}
var L = g;
const dr = /* @__PURE__ */ G(L), yr = /* @__PURE__ */ Q({
  __proto__: null,
  default: dr
}, [L]);
export {
  L as c,
  yr as i
};
