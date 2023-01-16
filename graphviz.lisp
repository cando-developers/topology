(in-package :foldamer)

(defun render-dag (dag stream)
  (format stream "digraph mygraph {~%")
  (maphash (lambda (key node)
             (declare (ignore key))
             (cond
               ((typep node 'cap-node)
                (format stream "  ~s~%" (symbol-name (name node))))
               ((typep node 'node)
                (format stream "  ~s~%" (symbol-name (name node))))
               ))
           (nodes dag))
  (loop for edge in (edges dag)
        for from-node = (from-node edge)
        for to-node = (to-node edge)
        for edge-type = (raw-name edge)
        do (format stream "~s -> ~s [label=~s]~%"
                   (symbol-name (name from-node))
                   (symbol-name (name to-node))
                   (symbol-name edge-type)))
  (format stream "}~%"))


(defun draw-dag (dag filename)
  (with-open-file (fout filename :direction :output)
    (render-dag dag fout)))


(defun node-name (node)
  (format nil "~a" (id node)))

(defun render-multiple-dags (dags stream)
  (format stream "digraph mygraph {~%")
  (loop for dag in dags
        for index from 0
        do (format stream "subgraph {~%")
        do (mapc (lambda (node)
                   (cond
                     ((eq node (root dag))
                      (format stream "  ~s[label=~s,shape=box,fillcolor=lightblue,style=filled]~%" (node-name node) (string (name node))))
                     ((typep node 'cap-node)
                      (format stream "  ~s[label=~s,shape=diamond]~%" (node-name node) (string (name node))))
                     ((typep node 'node)
                      (format stream "  ~s[label=~s]~%" (node-name node) (string (name node))))
                     ))
                 (nodes dag))
        do (loop for edge in (edges dag)
                 for from-node = (from-node edge)
                 for to-node = (to-node edge)
                 for edge-type = (raw-name edge)
                 do (format stream "  ~s -> ~s [label=~s]~%"
                            (node-name from-node)
                            (node-name to-node)
                            (symbol-name edge-type)))
        do (format stream "}~%")
        )
  (format stream "}~%"))

(defun draw-foldamer (foldamer filename)
  (with-open-file (fout filename :direction :output)
    (let ((dags (loop for training-oligomer-space in (foldamer::training-oligomer-spaces foldamer)
                      collect (foldamer::expression-dag training-oligomer-space))))
      (render-multiple-dags dags fout))))
