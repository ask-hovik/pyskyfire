pyskyfire.common.results
========================

.. py:module:: pyskyfire.common.results




Module Contents
---------------

.. py:class:: Results

   Bases: :py:obj:`collections.abc.MutableMapping`


   A mapping you can also dot-access, with save/load built in.


   .. py:attribute:: metadata


   .. py:method:: save(path)


   .. py:method:: load(path)
      :classmethod:



   .. py:method:: add(name, obj)

      Explicitly register something new.



   .. py:method:: list_items()


   .. py:method:: pop(key, default=__marker)

      D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
      If key is not found, d is returned if given, otherwise KeyError is raised.



   .. py:method:: popitem()

      D.popitem() -> (k, v), remove and return some (key, value) pair
      as a 2-tuple; but raise KeyError if D is empty.



   .. py:method:: clear()

      D.clear() -> None.  Remove all items from D.



   .. py:method:: update(other=(), /, **kwds)

      D.update([E, ]**F) -> None.  Update D from mapping/iterable E and F.
      If E present and has a .keys() method, does:     for k in E.keys(): D[k] = E[k]
      If E present and lacks .keys() method, does:     for (k, v) in E: D[k] = v
      In either case, this is followed by: for k, v in F.items(): D[k] = v



   .. py:method:: setdefault(key, default=None)

      D.setdefault(k[,d]) -> D.get(k,d), also set D[k]=d if k not in D



   .. py:method:: get(key, default=None)

      D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.



   .. py:method:: keys()

      D.keys() -> a set-like object providing a view on D's keys



   .. py:method:: items()

      D.items() -> a set-like object providing a view on D's items



   .. py:method:: values()

      D.values() -> an object providing a view on D's values



