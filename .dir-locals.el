
;; Emacs code style for files under this directory

;; For some reason, this doesn't work with cuda-mode.  So I instead use
;; c++-mode instead of cuda-mode for cu and cuh files.

((c++-mode . ((c-file-style . "Linux")
              (c-basic-offset . 4)
              (tab-width . 4)
              (indent-tabs-mode . t)
              (fill-column . 79))))
