import { Q as s } from "./index-8ee3baa4.js";
class h {
  constructor(e = {}) {
    if (!(e.maxSize && e.maxSize > 0))
      throw new TypeError("`maxSize` must be a number greater than 0");
    this.maxSize = e.maxSize, this.cache = /* @__PURE__ */ new Map(), this.oldCache = /* @__PURE__ */ new Map(), this._size = 0;
  }
  _set(e, t) {
    this.cache.set(e, t), this._size++, this._size >= this.maxSize && (this._size = 0, this.oldCache = this.cache, this.cache = /* @__PURE__ */ new Map());
  }
  get(e) {
    if (this.cache.has(e))
      return this.cache.get(e);
    if (this.oldCache.has(e)) {
      const t = this.oldCache.get(e);
      return this.oldCache.delete(e), this._set(e, t), t;
    }
  }
  set(e, t) {
    return this.cache.has(e) ? this.cache.set(e, t) : this._set(e, t), this;
  }
  has(e) {
    return this.cache.has(e) || this.oldCache.has(e);
  }
  peek(e) {
    if (this.cache.has(e))
      return this.cache.get(e);
    if (this.oldCache.has(e))
      return this.oldCache.get(e);
  }
  delete(e) {
    const t = this.cache.delete(e);
    return t && this._size--, this.oldCache.delete(e) || t;
  }
  clear() {
    this.cache.clear(), this.oldCache.clear(), this._size = 0;
  }
  *keys() {
    for (const [e] of this)
      yield e;
  }
  *values() {
    for (const [, e] of this)
      yield e;
  }
  *[Symbol.iterator]() {
    for (const e of this.cache)
      yield e;
    for (const e of this.oldCache) {
      const [t] = e;
      this.cache.has(t) || (yield e);
    }
  }
  get size() {
    let e = 0;
    for (const t of this.oldCache.keys())
      this.cache.has(t) || e++;
    return this._size + e;
  }
}
var i = h;
const r = /* @__PURE__ */ s(i);
export {
  r as Q
};
