import { Q as B, a1 as b } from "./index-8ee3baa4.js";
import { r as A } from "./core-2f49d2aa.js";
function O(p, h) {
  for (var f = 0; f < h.length; f++) {
    const a = h[f];
    if (typeof a != "string" && !Array.isArray(a)) {
      for (const c in a)
        if (c !== "default" && !(c in p)) {
          const i = Object.getOwnPropertyDescriptor(a, c);
          i && Object.defineProperty(p, c, i.get ? i : {
            enumerable: !0,
            get: () => a[c]
          });
        }
    }
  }
  return Object.freeze(Object.defineProperty(p, Symbol.toStringTag, { value: "Module" }));
}
var m = { exports: {} };
(function(p, h) {
  (function(f, a) {
    p.exports = a(A());
  })(b, function(f) {
    return function() {
      var a = f, c = a.lib, i = c.WordArray, x = a.enc;
      x.Base64 = {
        /**
         * Converts a word array to a Base64 string.
         *
         * @param {WordArray} wordArray The word array.
         *
         * @return {string} The Base64 string.
         *
         * @static
         *
         * @example
         *
         *     var base64String = CryptoJS.enc.Base64.stringify(wordArray);
         */
        stringify: function(t) {
          var s = t.words, n = t.sigBytes, o = this._map;
          t.clamp();
          for (var e = [], r = 0; r < n; r += 3)
            for (var v = s[r >>> 2] >>> 24 - r % 4 * 8 & 255, l = s[r + 1 >>> 2] >>> 24 - (r + 1) % 4 * 8 & 255, y = s[r + 2 >>> 2] >>> 24 - (r + 2) % 4 * 8 & 255, C = v << 16 | l << 8 | y, u = 0; u < 4 && r + u * 0.75 < n; u++)
              e.push(o.charAt(C >>> 6 * (3 - u) & 63));
          var g = o.charAt(64);
          if (g)
            for (; e.length % 4; )
              e.push(g);
          return e.join("");
        },
        /**
         * Converts a Base64 string to a word array.
         *
         * @param {string} base64Str The Base64 string.
         *
         * @return {WordArray} The word array.
         *
         * @static
         *
         * @example
         *
         *     var wordArray = CryptoJS.enc.Base64.parse(base64String);
         */
        parse: function(t) {
          var s = t.length, n = this._map, o = this._reverseMap;
          if (!o) {
            o = this._reverseMap = [];
            for (var e = 0; e < n.length; e++)
              o[n.charCodeAt(e)] = e;
          }
          var r = n.charAt(64);
          if (r) {
            var v = t.indexOf(r);
            v !== -1 && (s = v);
          }
          return _(t, s, o);
        },
        _map: "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="
      };
      function _(t, s, n) {
        for (var o = [], e = 0, r = 0; r < s; r++)
          if (r % 4) {
            var v = n[t.charCodeAt(r - 1)] << r % 4 * 2, l = n[t.charCodeAt(r)] >>> 6 - r % 4 * 2;
            o[e >>> 2] |= (v | l) << 24 - e % 4 * 8, e++;
          }
        return i.create(o, e);
      }
    }(), f.enc.Base64;
  });
})(m);
var d = m.exports;
const w = /* @__PURE__ */ B(d), D = /* @__PURE__ */ O({
  __proto__: null,
  default: w
}, [d]);
export {
  D as e
};
