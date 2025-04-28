import { a1 as k } from "./index-8ee3baa4.js";
var y = { exports: {} }, w;
function W() {
  return w || (w = 1, function(_, H) {
    (function(p, o) {
      _.exports = o();
    })(k, function() {
      var p = p || function(o, S) {
        var C = Object.create || function() {
          function t() {
          }
          return function(n) {
            var e;
            return t.prototype = n, e = new t(), t.prototype = null, e;
          };
        }(), d = {}, g = d.lib = {}, h = g.Base = function() {
          return {
            /**
             * Creates a new object that inherits from this object.
             *
             * @param {Object} overrides Properties to copy into the new object.
             *
             * @return {Object} The new object.
             *
             * @static
             *
             * @example
             *
             *     var MyType = CryptoJS.lib.Base.extend({
             *         field: 'value',
             *
             *         method: function () {
             *         }
             *     });
             */
            extend: function(t) {
              var n = C(this);
              return t && n.mixIn(t), (!n.hasOwnProperty("init") || this.init === n.init) && (n.init = function() {
                n.$super.init.apply(this, arguments);
              }), n.init.prototype = n, n.$super = this, n;
            },
            /**
             * Extends this object and runs the init method.
             * Arguments to create() will be passed to init().
             *
             * @return {Object} The new object.
             *
             * @static
             *
             * @example
             *
             *     var instance = MyType.create();
             */
            create: function() {
              var t = this.extend();
              return t.init.apply(t, arguments), t;
            },
            /**
             * Initializes a newly created object.
             * Override this method to add some logic when your objects are created.
             *
             * @example
             *
             *     var MyType = CryptoJS.lib.Base.extend({
             *         init: function () {
             *             // ...
             *         }
             *     });
             */
            init: function() {
            },
            /**
             * Copies properties into this object.
             *
             * @param {Object} properties The properties to mix in.
             *
             * @example
             *
             *     MyType.mixIn({
             *         field: 'value'
             *     });
             */
            mixIn: function(t) {
              for (var n in t)
                t.hasOwnProperty(n) && (this[n] = t[n]);
              t.hasOwnProperty("toString") && (this.toString = t.toString);
            },
            /**
             * Creates a copy of this object.
             *
             * @return {Object} The clone.
             *
             * @example
             *
             *     var clone = instance.clone();
             */
            clone: function() {
              return this.init.prototype.extend(this);
            }
          };
        }(), v = g.WordArray = h.extend({
          /**
           * Initializes a newly created word array.
           *
           * @param {Array} words (Optional) An array of 32-bit words.
           * @param {number} sigBytes (Optional) The number of significant bytes in the words.
           *
           * @example
           *
           *     var wordArray = CryptoJS.lib.WordArray.create();
           *     var wordArray = CryptoJS.lib.WordArray.create([0x00010203, 0x04050607]);
           *     var wordArray = CryptoJS.lib.WordArray.create([0x00010203, 0x04050607], 6);
           */
          init: function(t, n) {
            t = this.words = t || [], n != S ? this.sigBytes = n : this.sigBytes = t.length * 4;
          },
          /**
           * Converts this word array to a string.
           *
           * @param {Encoder} encoder (Optional) The encoding strategy to use. Default: CryptoJS.enc.Hex
           *
           * @return {string} The stringified word array.
           *
           * @example
           *
           *     var string = wordArray + '';
           *     var string = wordArray.toString();
           *     var string = wordArray.toString(CryptoJS.enc.Utf8);
           */
          toString: function(t) {
            return (t || b).stringify(this);
          },
          /**
           * Concatenates a word array to this word array.
           *
           * @param {WordArray} wordArray The word array to append.
           *
           * @return {WordArray} This word array.
           *
           * @example
           *
           *     wordArray1.concat(wordArray2);
           */
          concat: function(t) {
            var n = this.words, e = t.words, r = this.sigBytes, s = t.sigBytes;
            if (this.clamp(), r % 4)
              for (var i = 0; i < s; i++) {
                var f = e[i >>> 2] >>> 24 - i % 4 * 8 & 255;
                n[r + i >>> 2] |= f << 24 - (r + i) % 4 * 8;
              }
            else
              for (var i = 0; i < s; i += 4)
                n[r + i >>> 2] = e[i >>> 2];
            return this.sigBytes += s, this;
          },
          /**
           * Removes insignificant bits.
           *
           * @example
           *
           *     wordArray.clamp();
           */
          clamp: function() {
            var t = this.words, n = this.sigBytes;
            t[n >>> 2] &= 4294967295 << 32 - n % 4 * 8, t.length = o.ceil(n / 4);
          },
          /**
           * Creates a copy of this word array.
           *
           * @return {WordArray} The clone.
           *
           * @example
           *
           *     var clone = wordArray.clone();
           */
          clone: function() {
            var t = h.clone.call(this);
            return t.words = this.words.slice(0), t;
          },
          /**
           * Creates a word array filled with random bytes.
           *
           * @param {number} nBytes The number of random bytes to generate.
           *
           * @return {WordArray} The random word array.
           *
           * @static
           *
           * @example
           *
           *     var wordArray = CryptoJS.lib.WordArray.random(16);
           */
          random: function(t) {
            for (var n = [], e = function(a) {
              var a = a, c = 987654321, u = 4294967295;
              return function() {
                c = 36969 * (c & 65535) + (c >> 16) & u, a = 18e3 * (a & 65535) + (a >> 16) & u;
                var l = (c << 16) + a & u;
                return l /= 4294967296, l += 0.5, l * (o.random() > 0.5 ? 1 : -1);
              };
            }, r = 0, s; r < t; r += 4) {
              var i = e((s || o.random()) * 4294967296);
              s = i() * 987654071, n.push(i() * 4294967296 | 0);
            }
            return new v.init(n, t);
          }
        }), x = d.enc = {}, b = x.Hex = {
          /**
           * Converts a word array to a hex string.
           *
           * @param {WordArray} wordArray The word array.
           *
           * @return {string} The hex string.
           *
           * @static
           *
           * @example
           *
           *     var hexString = CryptoJS.enc.Hex.stringify(wordArray);
           */
          stringify: function(t) {
            for (var n = t.words, e = t.sigBytes, r = [], s = 0; s < e; s++) {
              var i = n[s >>> 2] >>> 24 - s % 4 * 8 & 255;
              r.push((i >>> 4).toString(16)), r.push((i & 15).toString(16));
            }
            return r.join("");
          },
          /**
           * Converts a hex string to a word array.
           *
           * @param {string} hexStr The hex string.
           *
           * @return {WordArray} The word array.
           *
           * @static
           *
           * @example
           *
           *     var wordArray = CryptoJS.enc.Hex.parse(hexString);
           */
          parse: function(t) {
            for (var n = t.length, e = [], r = 0; r < n; r += 2)
              e[r >>> 3] |= parseInt(t.substr(r, 2), 16) << 24 - r % 8 * 4;
            return new v.init(e, n / 2);
          }
        }, B = x.Latin1 = {
          /**
           * Converts a word array to a Latin1 string.
           *
           * @param {WordArray} wordArray The word array.
           *
           * @return {string} The Latin1 string.
           *
           * @static
           *
           * @example
           *
           *     var latin1String = CryptoJS.enc.Latin1.stringify(wordArray);
           */
          stringify: function(t) {
            for (var n = t.words, e = t.sigBytes, r = [], s = 0; s < e; s++) {
              var i = n[s >>> 2] >>> 24 - s % 4 * 8 & 255;
              r.push(String.fromCharCode(i));
            }
            return r.join("");
          },
          /**
           * Converts a Latin1 string to a word array.
           *
           * @param {string} latin1Str The Latin1 string.
           *
           * @return {WordArray} The word array.
           *
           * @static
           *
           * @example
           *
           *     var wordArray = CryptoJS.enc.Latin1.parse(latin1String);
           */
          parse: function(t) {
            for (var n = t.length, e = [], r = 0; r < n; r++)
              e[r >>> 2] |= (t.charCodeAt(r) & 255) << 24 - r % 4 * 8;
            return new v.init(e, n);
          }
        }, z = x.Utf8 = {
          /**
           * Converts a word array to a UTF-8 string.
           *
           * @param {WordArray} wordArray The word array.
           *
           * @return {string} The UTF-8 string.
           *
           * @static
           *
           * @example
           *
           *     var utf8String = CryptoJS.enc.Utf8.stringify(wordArray);
           */
          stringify: function(t) {
            try {
              return decodeURIComponent(escape(B.stringify(t)));
            } catch {
              throw new Error("Malformed UTF-8 data");
            }
          },
          /**
           * Converts a UTF-8 string to a word array.
           *
           * @param {string} utf8Str The UTF-8 string.
           *
           * @return {WordArray} The word array.
           *
           * @static
           *
           * @example
           *
           *     var wordArray = CryptoJS.enc.Utf8.parse(utf8String);
           */
          parse: function(t) {
            return B.parse(unescape(encodeURIComponent(t)));
          }
        }, m = g.BufferedBlockAlgorithm = h.extend({
          /**
           * Resets this block algorithm's data buffer to its initial state.
           *
           * @example
           *
           *     bufferedBlockAlgorithm.reset();
           */
          reset: function() {
            this._data = new v.init(), this._nDataBytes = 0;
          },
          /**
           * Adds new data to this block algorithm's buffer.
           *
           * @param {WordArray|string} data The data to append. Strings are converted to a WordArray using UTF-8.
           *
           * @example
           *
           *     bufferedBlockAlgorithm._append('data');
           *     bufferedBlockAlgorithm._append(wordArray);
           */
          _append: function(t) {
            typeof t == "string" && (t = z.parse(t)), this._data.concat(t), this._nDataBytes += t.sigBytes;
          },
          /**
           * Processes available data blocks.
           *
           * This method invokes _doProcessBlock(offset), which must be implemented by a concrete subtype.
           *
           * @param {boolean} doFlush Whether all blocks and partial blocks should be processed.
           *
           * @return {WordArray} The processed data.
           *
           * @example
           *
           *     var processedData = bufferedBlockAlgorithm._process();
           *     var processedData = bufferedBlockAlgorithm._process(!!'flush');
           */
          _process: function(t) {
            var n = this._data, e = n.words, r = n.sigBytes, s = this.blockSize, i = s * 4, f = r / i;
            t ? f = o.ceil(f) : f = o.max((f | 0) - this._minBufferSize, 0);
            var a = f * s, c = o.min(a * 4, r);
            if (a) {
              for (var u = 0; u < a; u += s)
                this._doProcessBlock(e, u);
              var l = e.splice(0, a);
              n.sigBytes -= c;
            }
            return new v.init(l, c);
          },
          /**
           * Creates a copy of this object.
           *
           * @return {Object} The clone.
           *
           * @example
           *
           *     var clone = bufferedBlockAlgorithm.clone();
           */
          clone: function() {
            var t = h.clone.call(this);
            return t._data = this._data.clone(), t;
          },
          _minBufferSize: 0
        });
        g.Hasher = m.extend({
          /**
           * Configuration options.
           */
          cfg: h.extend(),
          /**
           * Initializes a newly created hasher.
           *
           * @param {Object} cfg (Optional) The configuration options to use for this hash computation.
           *
           * @example
           *
           *     var hasher = CryptoJS.algo.SHA256.create();
           */
          init: function(t) {
            this.cfg = this.cfg.extend(t), this.reset();
          },
          /**
           * Resets this hasher to its initial state.
           *
           * @example
           *
           *     hasher.reset();
           */
          reset: function() {
            m.reset.call(this), this._doReset();
          },
          /**
           * Updates this hasher with a message.
           *
           * @param {WordArray|string} messageUpdate The message to append.
           *
           * @return {Hasher} This hasher.
           *
           * @example
           *
           *     hasher.update('message');
           *     hasher.update(wordArray);
           */
          update: function(t) {
            return this._append(t), this._process(), this;
          },
          /**
           * Finalizes the hash computation.
           * Note that the finalize operation is effectively a destructive, read-once operation.
           *
           * @param {WordArray|string} messageUpdate (Optional) A final message update.
           *
           * @return {WordArray} The hash.
           *
           * @example
           *
           *     var hash = hasher.finalize();
           *     var hash = hasher.finalize('message');
           *     var hash = hasher.finalize(wordArray);
           */
          finalize: function(t) {
            t && this._append(t);
            var n = this._doFinalize();
            return n;
          },
          blockSize: 16,
          /**
           * Creates a shortcut function to a hasher's object interface.
           *
           * @param {Hasher} hasher The hasher to create a helper for.
           *
           * @return {Function} The shortcut function.
           *
           * @static
           *
           * @example
           *
           *     var SHA256 = CryptoJS.lib.Hasher._createHelper(CryptoJS.algo.SHA256);
           */
          _createHelper: function(t) {
            return function(n, e) {
              return new t.init(e).finalize(n);
            };
          },
          /**
           * Creates a shortcut function to the HMAC's object interface.
           *
           * @param {Hasher} hasher The hasher to use in this HMAC helper.
           *
           * @return {Function} The shortcut function.
           *
           * @static
           *
           * @example
           *
           *     var HmacSHA256 = CryptoJS.lib.Hasher._createHmacHelper(CryptoJS.algo.SHA256);
           */
          _createHmacHelper: function(t) {
            return function(n, e) {
              return new F.HMAC.init(t, e).finalize(n);
            };
          }
        });
        var F = d.algo = {};
        return d;
      }(Math);
      return p;
    });
  }(y)), y.exports;
}
export {
  W as r
};
