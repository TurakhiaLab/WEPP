import { Q as N, a1 as T } from "./index-8ee3baa4.js";
import { r as $ } from "./core-2f49d2aa.js";
function G(p, A) {
  for (var c = 0; c < A.length; c++) {
    const t = A[c];
    if (typeof t != "string" && !Array.isArray(t)) {
      for (const n in t)
        if (n !== "default" && !(n in p)) {
          const h = Object.getOwnPropertyDescriptor(t, n);
          h && Object.defineProperty(p, n, h.get ? h : {
            enumerable: !0,
            get: () => t[n]
          });
        }
    }
  }
  return Object.freeze(Object.defineProperty(p, Symbol.toStringTag, { value: "Module" }));
}
var S = { exports: {} };
(function(p, A) {
  (function(c, t) {
    p.exports = t($());
  })(T, function(c) {
    return function(t) {
      var n = c, h = n.lib, C = h.WordArray, H = h.Hasher, W = n.algo, w = [], B = [];
      (function() {
        function i(o) {
          for (var l = t.sqrt(o), u = 2; u <= l; u++)
            if (!(o % u))
              return !1;
          return !0;
        }
        function v(o) {
          return (o - (o | 0)) * 4294967296 | 0;
        }
        for (var r = 2, e = 0; e < 64; )
          i(r) && (e < 8 && (w[e] = v(t.pow(r, 1 / 2))), B[e] = v(t.pow(r, 1 / 3)), e++), r++;
      })();
      var f = [], j = W.SHA256 = H.extend({
        _doReset: function() {
          this._hash = new C.init(w.slice(0));
        },
        _doProcessBlock: function(i, v) {
          for (var r = this._hash.words, e = r[0], o = r[1], l = r[2], u = r[3], s = r[4], y = r[5], b = r[6], g = r[7], a = 0; a < 64; a++) {
            if (a < 16)
              f[a] = i[v + a] | 0;
            else {
              var _ = f[a - 15], q = (_ << 25 | _ >>> 7) ^ (_ << 14 | _ >>> 18) ^ _ >>> 3, d = f[a - 2], x = (d << 15 | d >>> 17) ^ (d << 13 | d >>> 19) ^ d >>> 10;
              f[a] = q + f[a - 7] + x + f[a - 16];
            }
            var D = s & y ^ ~s & b, F = e & o ^ e & l ^ o & l, k = (e << 30 | e >>> 2) ^ (e << 19 | e >>> 13) ^ (e << 10 | e >>> 22), z = (s << 26 | s >>> 6) ^ (s << 21 | s >>> 11) ^ (s << 7 | s >>> 25), P = g + z + D + B[a] + f[a], E = k + F;
            g = b, b = y, y = s, s = u + P | 0, u = l, l = o, o = e, e = P + E | 0;
          }
          r[0] = r[0] + e | 0, r[1] = r[1] + o | 0, r[2] = r[2] + l | 0, r[3] = r[3] + u | 0, r[4] = r[4] + s | 0, r[5] = r[5] + y | 0, r[6] = r[6] + b | 0, r[7] = r[7] + g | 0;
        },
        _doFinalize: function() {
          var i = this._data, v = i.words, r = this._nDataBytes * 8, e = i.sigBytes * 8;
          return v[e >>> 5] |= 128 << 24 - e % 32, v[(e + 64 >>> 9 << 4) + 14] = t.floor(r / 4294967296), v[(e + 64 >>> 9 << 4) + 15] = r, i.sigBytes = v.length * 4, this._process(), this._hash;
        },
        clone: function() {
          var i = H.clone.call(this);
          return i._hash = this._hash.clone(), i;
        }
      });
      n.SHA256 = H._createHelper(j), n.HmacSHA256 = H._createHmacHelper(j);
    }(Math), c.SHA256;
  });
})(S);
var O = S.exports;
const K = /* @__PURE__ */ N(O), R = /* @__PURE__ */ G({
  __proto__: null,
  default: K
}, [O]);
export {
  R as s
};
