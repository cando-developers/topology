(in-package :foldamer)


(defun done-data (trainer-context)
  (multiple-value-bind (input-file done-file)
    (foldamer:calculate-files trainer-context)
    (declare (ignore input-file))
    (if (null (probe-file done-file))
	(values nil)
      (with-open-file (done-fin done-file :direction :input)
        (read done-fin)))))

(defun needs-work (trainer-context steps)
  "Return (values needs-work hits) 
    If this trainer-context has not run enough steps then needs-work is T.
    Also return the number of hits that have been found."
  (let* ((done-data (done-data trainer-context))
	 (done-steps (getf done-data :finished-steps))
	 (hits (getf done-data :hits)))
    (or (null done-data) (or (< done-steps steps) (= hits 0)))))

(defun max-finished-steps (spiros)
  (let ((max-finished-steps 0)
        (valid-trainer-contexts (foldamer:valid-trainer-contexts spiros)))
    (loop for context in valid-trainer-contexts
          do (multiple-value-bind (input-file done-file)
                 (foldamer:calculate-files context)
               (declare (ignore input-file))
               #+(or)(format t "Reading ~a~%" done-file)
               (when (probe-file done-file)
                 (with-open-file (fin done-file)
                   (let ((info (read fin)))
                     (setf max-finished-steps (max max-finished-steps (getf info :finished-steps)))
                     )))))
    max-finished-steps))


(defun needs-extract-conformations (path spiros)
  (let ((valid-trainer-contexts (foldamer:valid-trainer-contexts spiros))
        (conformations-pathname (merge-pathnames *conformations-filename* path)))
    (if (null (probe-file conformations-pathname))
        t ;; Need to update conformations if file doesn't exist
        (let ((latest-done-file-write-date 0)
              (conformations-write-date (file-write-date conformations-pathname)))
          (loop for context in valid-trainer-contexts
                do (multiple-value-bind (input-file done-file)
                       (foldamer:calculate-files context)
                     (declare (ignore input-file))
                     #+(or)(format t "Reading ~a~%" done-file)
                     (when (probe-file done-file)
                       (setf latest-done-file-write-date (max latest-done-file-write-date (file-write-date done-file))))))
          (< conformations-write-date latest-done-file-write-date)))))

(defun foldamer-update-trainers (spiros &key remove-unused force-save)
  (foldamer:maybe-remove-unused-trainers spiros #P"./" :doit remove-unused)
  (foldamer:generate-training-oligomers spiros #P"./" :force-save force-save))

(defun foldamer-status (&key spiros)
  (let* ((trainer-contexts (foldamer:valid-trainer-contexts spiros))
	 (trainer-done-data (make-hash-table :test 'equal))
	 (max-finished-steps (max-finished-steps spiros)))
    (loop for trainer-context in trainer-contexts
	  for done-data = (done-data trainer-context)
	  do (setf (gethash trainer-context trainer-done-data) done-data))
    (loop for trainer-context in trainer-contexts
	  for done-data = (gethash trainer-context trainer-done-data)
	  for hits = (getf done-data :hits)
	  for finished-steps = (getf done-data :finished-steps)
	  when (and hits (= hits 0))
            do (format t "~a :finished-steps ~a :hits ~a~%" trainer-context finished-steps hits))
    (let ((number-of-trainers-with-max-finished-steps 
	   (loop for trainer-context in trainer-contexts
		 for done-data = (gethash trainer-context trainer-done-data)
		 for finished-steps = (getf done-data :finished-steps)
		 when (and done-data (= finished-steps max-finished-steps))
		   count 1)))
      (format t "Number of trainers with ~a finished steps = ~a~%" 
	      max-finished-steps
	      number-of-trainers-with-max-finished-steps))
    (format t "Number of trainers ~a~%" (length trainer-contexts))
    ))
			      
(defun foldamer-setup (spiros
                         &key
			 (nodes 12)
                         (clasp "$CLASP -l code.lisp")
                         (script "jobs.dat")
                         (finished-steps 2 finished-steps-p)
                         (advance-steps 10 advance-steps-p))
  (foldamer-update-trainers spiros :remove-unused t)
  (let* ((max-finished-steps (max-finished-steps spiros))
         (steps (cond
                  ((and (null advance-steps-p) finished-steps-p) finished-steps)
                  ((and (null finished-steps-p) advance-steps-p) (+ max-finished-steps advance-steps))
                  ((and (null finished-steps-p) (null advance-steps-p))
                   max-finished-steps)
                  (t (error "You must provide only one option advance-steps or finished-steps or neither"))))
         (trainer-contexts (foldamer:valid-trainer-contexts spiros))
         (node-total nodes)
	 (jobs-per-node (ceiling (length trainer-contexts) node-total))
         (node-work (make-hash-table :test 'eql))
         (trainer-count 0)
         (trainers-that-need-work (loop for trainer-context in trainer-contexts
                                        when (needs-work trainer-context steps)
                                          collect trainer-context)))
    
    (format t "Maximum finished-steps ~a~%" max-finished-steps)
    (format t "New steps              ~a~%" steps)
    (loop for trainer-file in trainers-that-need-work
          for index from 0
          for node-index = (mod index node-total)
          for trainer-context = (pathname-name (pathname trainer-file))
          do (progn
               (push trainer-context (gethash node-index node-work nil))
               (incf trainer-count)
               (when (gethash node-index node-work)
                 #+(or)(format t "~a ~a~%" node-index trainer-context))))
    (with-open-file (fout script :direction :output :if-exists :supersede)
      (loop for node-index below (min jobs-per-node (hash-table-count node-work))
            do (format fout "~a -e \"(foldamer:foldamer-run-node ~a :steps ~a)\"~%"
                       clasp
                       node-index
                       steps)))
    (cando:save-cando node-work "node-work.cando")
    (format t "~a trainers need to be trained~%" (if (= trainer-count 0) "No" trainer-count))
    (when (= trainer-count 0)
      (if (needs-extract-conformations #P"./" spiros)
          (format t "Needs foldamer-extract-conformations~%")
          (format t "Ready to build spiroligomers~%"))))
  (sys:quit))


(defun foldamer-extract-conformations (&optional (path #P"./"))
  (let ((foldamer-conformations-map (foldamer::assemble-fragment-conformations-map path)))
    (format t "Saving ~a~%" *conformations-filename*)
    (cando:save-cando foldamer-conformations-map *conformations-filename*)))


(defun foldamer-run-node (node-index &key (steps 2) testing (parallel t) (no-fail t))
  (format t "Loading force-field~%")
  (foldamer:prepare-to-build-trainer :smirnoff "./force-field.offxml")
  (let* ((all-node-work (cando:load-cando "node-work.cando"))
         (node-trainers (gethash node-index all-node-work))
         (foldamer-dat-pathname (merge-pathnames #P"./foldamer.dat"))
         (foldamer (cando:load-cando foldamer-dat-pathname)))
    (flet ((one-trainer (trainer-context)
             (let ((input-file (merge-pathnames (make-pathname :directory '(:relative "data")
                                                               :name trainer-context
                                                               :type "input"))))
               (format t "About to build-trainer: ~a ~a~%" trainer-context (probe-file input-file))
               (if no-fail
                   (handler-bind
                       ((error (lambda (err)
                                 (format t "~s~%" err)
                                 (let ((clasp-debug:*frame-filters* nil))
                                   (clasp-debug:print-backtrace)))))
                     (foldamer:build-trainer foldamer trainer-context :load-pathname input-file :steps steps))
                   (foldamer:build-trainer foldamer trainer-context :load-pathname input-file :steps steps)))))
      (if parallel
          (lparallel:pmapc #'one-trainer node-trainers)
          (mapc #'one-trainer node-trainers))))
  (unless testing (sys:quit)))

