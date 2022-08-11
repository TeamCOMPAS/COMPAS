.. things that are needed globally ...  JR


.. newline - inserts '<br />' in the html

.. |br| raw:: html

   <br />


.. nbsp - insertt a non-breaking space in the html

.. |_| unicode:: 0xA0 
   :trim:


.. rst doesn't seem to support nesting decorations (e.g. can't do '**text**' to get 'text' in bold inside quote marks
.. this implements bold text and italic text via CSS (along with the corresponding entry in the CSS file) which can then
.. be wrapped in quote marks (or whatever else you want to wrap it in) e.g. some text ':boldtext:`bold_text`'
.. (I wanted to just use role :bold: and :italic:, but I think at least :bold: is defined elsewhere in readthedocs... 
.. so I made :bolditalictext: consistent)

.. role:: boldtext
  :class: boldtext

.. role:: italictext
  :class: italictext


.. rst doesn't support bold italics - this implements it via CSS (along with the corresponding entry in the CSS file)

.. role:: bolditalictext
  :class: bolditalictext

