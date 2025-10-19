pyskyfire.common.results.Results
================================

.. py:class:: pyskyfire.common.results.Results

   Bases: :py:obj:`collections.abc.MutableMapping`

   .. autoapi-inheritance-diagram:: pyskyfire.common.results.Results
      :parts: 1
      :private-bases:


   
   Dictionary-like container with dot access and persistence helpers.

   Behaves as a normal mutable mapping while allowing attribute-style
   access (``res.foo`` ≡ ``res["foo"]``) and built-in save/load using
   :mod:`cloudpickle`.


   :Attributes:

       **_store** : :class:`python:dict`\[:class:`python:str`, :obj:`Any`]
           Internal data store for all key/value pairs.

       **metadata** : :class:`python:dict`\[:class:`python:str`, :obj:`Any`]
           Optional auxiliary metadata (e.g., creation date, simulation info).









   .. seealso::

       
       :obj:`dict`
           Built-in Python mapping type.
       :obj:`cloudpickle`
           Used for serialization in :meth:`save` and :meth:`load`.
       
       
   .. rubric:: Notes

   Dot assignment automatically writes into :attr:`_store`. Internal
   attributes such as ``_store`` and ``metadata`` are protected from
   accidental overwrite.



   ..
       !! processed by numpydoc !!

   .. py:method:: __delitem__(key)

      
      Delete ``key`` from the store.
















      ..
          !! processed by numpydoc !!


   .. py:method:: __getattr__(name)

      
      Return an attribute or stored item.


      :Parameters:

          **name** : :class:`python:str`
              Attribute or key name.



      :Returns:

          :obj:`Any`
              The stored object corresponding to ``name``.




      :Raises:

          :obj:`AttributeError`
              If ``name`` is not found in :attr:`_store` or the instance.







      ..
          !! processed by numpydoc !!


   .. py:method:: __getitem__(key)

      
      Return the value associated with ``key``.
















      ..
          !! processed by numpydoc !!


   .. py:method:: __iter__()

      
      Iterate over stored keys.
















      ..
          !! processed by numpydoc !!


   .. py:method:: __len__()

      
      Return the number of stored items.
















      ..
          !! processed by numpydoc !!


   .. py:method:: __setattr__(name, value)

      
      Assign attributes to the internal store by default.


      :Parameters:

          **name** : :class:`python:str`
              Attribute name.

          **value** : :obj:`Any`
              Value to store or assign.











      .. rubric:: Notes

      ``_store`` and ``metadata`` are written directly to ``__dict__``;
      all other names go into :attr:`_store`.



      ..
          !! processed by numpydoc !!


   .. py:method:: __setitem__(key, value)

      
      Store ``value`` under ``key``.
















      ..
          !! processed by numpydoc !!


   .. py:method:: add(name, obj)

      
      Explicitly register a new object.


      :Parameters:

          **name** : :class:`python:str`
              Key under which to store the object.

          **obj** : :obj:`Any`
              The object to add.














      ..
          !! processed by numpydoc !!


   .. py:method:: list_items()

      
      List all stored keys.





      :Returns:

          :class:`python:list`\[:class:`python:str`]
              Keys currently stored in the container.











      ..
          !! processed by numpydoc !!


   .. py:method:: load(path)
      :classmethod:


      
      Load a previously saved :class:`Results` object.


      :Parameters:

          **path** : :class:`python:str` | :obj:`PathLike`
              Path to the pickled file.



      :Returns:

          :obj:`Results`
              The deserialized object instance.








      .. rubric:: Notes

      Uses :func:`cloudpickle.load` for unpickling.



      ..
          !! processed by numpydoc !!


   .. py:method:: save(path)

      
      Serialize this :class:`Results` object to disk.


      :Parameters:

          **path** : :class:`python:str` | :obj:`PathLike`
              File path to write the pickled object.











      .. rubric:: Notes

      Uses :func:`cloudpickle.dump` for robust serialization of arbitrary
      Python objects.



      ..
          !! processed by numpydoc !!

